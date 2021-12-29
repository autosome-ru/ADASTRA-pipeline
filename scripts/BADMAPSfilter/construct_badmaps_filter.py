import json
from multiprocessing import Process, Pool

import pandas as pd
from pandas.errors import EmptyDataError
import os
import numpy as np
from scipy.stats import levene
import statsmodels.stats.multitest

from tqdm import tqdm

from scripts.HELPERS.helpers import pack, segmentation_states, get_babachi_models_list, get_states_from_model_name
from scripts.HELPERS.paths import get_excluded_badmaps_list_path, get_correlation_file_path, get_correlation_path, \
    get_release_stats_path

from scipy import stats as st
from joblib import Parallel, delayed


def find_test_datasets(cor_stats_path):
    df = pd.read_table(cor_stats_path)
    return df


def get_path(row, model, remake=False):
    return os.path.join(get_correlation_path(), '{}_tables{}/'.format(model, '_filtered' if remake else '')
                              + row['#cell_line'] + '@' + row['cells'] + '.tsv')


def get_data_from_cor_row(row, model, remake=False):
    group = make_dataset(row['#cell_line'], row['cells'])
    path = get_path(row, model, remake=remake)
    try:
        tmp_df = pd.read_table(path, comment='#', header=None)
    except EmptyDataError:
        return False, None, None
    return True, group, tmp_df


def open_dfs(df, model, concat=True, remake=False):
    res_dfs = []
    for index, row in df.iterrows():
        ok, group, tmp_df = get_data_from_cor_row(row, model, remake=remake)
        if not ok:
            continue
        states = get_states_from_model_name(model)
        tmp_df.columns = ['chr', 'pos', 'ref', 'alt', 'BAD'] + ['Q{:.2f}'.format(b) for b in
                                                               states] + ['snps', 'cov',
                                                                                       'dataset',
                                                                                       'seg_id',
                                                                                       'p_value']
        if not concat:
            tmp_df['cov'] = tmp_df['ref'] + tmp_df['alt']
            tmp_df['max'] = tmp_df[['ref', 'alt']].max(axis=1)
            tmp_df['min'] = tmp_df[['ref', 'alt']].min(axis=1)
            # tmp_df['es'] = -np.log2(tmp_df['max'] / tmp_df['min'] / tmp_df['BAD'])

            stats = {
                BAD: collect_stats_df(tmp_df, BAD).to_dict()
                for BAD in states
            }

            if res_dfs is None:
                res_dfs = [(group, stats)]
            else:
                res_dfs.append((group, stats))
        else:
            tmp_df['group'] = group
            tmp_df['cell_line'] = row['#cell_line']
            if res_dfs is None:
                res_dfs = [tmp_df]
            else:
                res_dfs.append(tmp_df)
    if concat:
        # FIXME collect stats here too
        res_df = pd.concat(res_dfs)
        res_df['cov'] = res_df['ref'] + res_df['alt']
        res_df['max'] = res_df[['ref', 'alt']].max(axis=1)
        res_df['min'] = res_df[['ref', 'alt']].min(axis=1)
        # res_df['es'] = -np.log2(res_df['max'] / res_df['min'] / res_df['BAD'])
        return res_df
    else:
        return res_dfs


def make_dataset(cell_line, lab):
    return '{}@{}'.format(cell_line, lab)


def collect_stats_df(df, BAD):
    out_t = df[df['BAD'] == BAD].groupby(['ref', 'alt']).size().reset_index(name='counts')
    out_t.fillna(0, inplace=True)
    out_t.columns = ['ref', 'alt', 'counts']
    return out_t


def rmsea_gof(stat, df, norm):
    if norm <= 1:
        return 0
    else:
        # if max(stat - df, 0) / (df * (norm - 1)) < 0:
        #     print(stat, df)
        score = np.sqrt(max(stat - df, 0) / (df * (norm - 1)))
    return score


def calculate_gof(counts_array, expected, norm, number_of_params):
    observed = counts_array.copy()

    idxs = (observed != 0) & (expected != 0)
    if idxs.sum() <= number_of_params + 1:
        return 0
    df = idxs.sum() - 1 - number_of_params
    stat = np.sum(observed[idxs] * (np.log(observed[idxs]) - np.log(expected[idxs]))) * 2
    return rmsea_gof(stat, df, norm)


def make_binom_density(cov, BAD, allele_tr):
    b1 = st.binom(cov, 1 / (BAD + 1))
    b2 = st.binom(cov, BAD / (BAD + 1))
    pdf1 = b1.pmf
    pdf2 = b2.pmf
    cdf = lambda x: 0.5 * (b1.cdf(x) + b2.cdf(x))
    sf = lambda x: 0.5 * (b1.sf(x) + b2.sf(x))
    norm = cdf(cov - allele_tr) + sf(allele_tr - 1) - 1

    res = np.zeros(cov + 1, dtype=np.float_)
    for i in range(allele_tr, cov - allele_tr + 1):
        res[i] = 0.5 * (pdf1(i) + pdf2(i))
    return res / norm


def init_process_for_mode(args):
    mode, cors = args
    test_dfs = open_dfs(cors, mode, remake=False, concat=False)
    print('Test concatenated {}'.format(mode))
    states = get_states_from_model_name(mode)

    with open(get_excluded_badmaps_list_path(model=mode, remake=False), 'w') as out:
        out.write(
            pack(['#cell_line', 'sample', 'correlation'] + ['RMSEA_GOF_BAD{:.2f}'.format(state) for state in states] +
                 ['NUM_TESTED_SNPS_BAD{:.2f}'.format(state) for state in states]))

    return test_dfs


def process_for_dataset(mode, dataset, cors, min_cov, max_cov):
    states = get_states_from_model_name(mode)
    try:
        cell_line, lab = dataset.split('@')
    except Exception:
        print(dataset)
        raise
    cor = cors[(cors['#cell_line'] == cell_line) & (cors['cells'] == lab)][
        'cor_by_snp_{}'.format(mode)].tolist()
    assert len(cor) == 1
    cor = cor[0]

    stats_dir = os.path.join(get_release_stats_path(), 'filter_stats')
    file_name = '{}_{}_stats.json'.format(dataset, mode)
    print(file_name)
    with open(os.path.join(stats_dir, file_name)) as file:
        stats_for_bads_dict = json.load(file)

    gofs = {}
    tested_snps = {}
    for BAD in states:
        stats = pd.DataFrame(stats_for_bads_dict[str(BAD)], dtype=np.int_)
        counts_arrays = []
        expected_arrays = []
        for cov in range(min_cov, max_cov + 1):
            counts_array = np.zeros(cov + 1, dtype=np.int_)
            for ref in range(5, cov - 5 + 1):
                counts = stats[(stats['ref'] == ref) & (stats['alt'] == cov - ref)]['counts'].to_list()
                if len(counts) == 1:
                    counts_array[ref] = counts[0]
                elif len(counts) > 1:
                    print(counts)
            counts_arrays.append(counts_array)
            expected_arrays.append(make_binom_density(cov, BAD, 5) * counts_array.sum())
        observed = np.concatenate(counts_arrays)
        expected = np.concatenate(expected_arrays)
        norm = observed.sum()
        tested_snps[BAD] = norm
        gofs[BAD] = calculate_gof(observed, expected, norm, 1)

    with open(get_excluded_badmaps_list_path(model=mode, remake=False), 'a') as out:
        out.write(pack([cell_line, lab, cor] +
                       [gofs[BAD] for BAD in states] +
                       [tested_snps[BAD] for BAD in states]))


def main(min_cov, max_cov, n_jobs, collect_stats=True):
    modes = get_babachi_models_list(remake=False)
    correlation_file_path = get_correlation_file_path(remake=False)
    cors = pd.read_table(correlation_file_path)
    stats_dir = os.path.join(get_release_stats_path(), 'filter_stats')

    print('Preprocessing started')
    if collect_stats:
        pool = Pool(processes=5)
        test_dfs_lists = pool.map(
            init_process_for_mode, ((mode, cors) for mode in modes)
        )

        if not os.path.isdir(stats_dir):
            os.mkdir(stats_dir)
        for mode, test_dfs in zip(modes, test_dfs_lists):
            for dataset, stats in test_dfs:
                file_name = '{}_{}_stats.json'.format(dataset, mode)
                with open(os.path.join(stats_dir, file_name), 'w') as file:
                    json.dump(stats, file, indent=2)
        datasets_lists = [
            [dataset for dataset, _ in test_dfs]
            for test_dfs in test_dfs_lists
        ]

    else:
        datasets_lists = []
        for mode in modes:
            list_for_mode = []
            for index, row in tqdm(cors.iterrows()):
                ok, dataset, _ = get_data_from_cor_row(row, mode)
                if ok:
                    list_for_mode.append(dataset)
            datasets_lists.append(list_for_mode)

    def dataset_process(mode, dataset):
        process_for_dataset(mode, dataset, cors, min_cov, max_cov)

    print('Preprocessing finished')
    Parallel(n_jobs=n_jobs, verbose=10)(delayed(dataset_process)(mode, dataset)
                       for mode, datasets in zip(modes, datasets_lists)
                       for dataset in datasets)
