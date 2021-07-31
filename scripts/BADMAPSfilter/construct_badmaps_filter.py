import json

import pandas as pd
from pandas.errors import EmptyDataError
import os
import numpy as np
from scipy.stats import percentileofscore, levene
import statsmodels.stats.multitest

from scripts.HELPERS.helpers import pack
from scripts.HELPERS.paths import get_excluded_badmaps_list_path, get_correlation_file_path

big_cell_lines = [
    'K562__myelogenous_leukemia_',
    'MCF7__Invasive_ductal_breast_carcinoma_',
    'A549__lung_carcinoma_',
    '22RV1__prostate_carcinoma_',
    'HCT-116__colon_carcinoma_',
]


def find_ref_datasets(cor_stats_path):
    df = pd.read_table(cor_stats_path)
    df = df[df['cor_by_snp_CAIC'] > 0.75]
    return df


def find_test_datasets(cor_stats_path):
    df = pd.read_table(cor_stats_path)
    return df


def get_path(row, remake=False):
    return os.path.expanduser('~/DataChipInt/BADmaps/Correlation/CAIC_tables{}/'.format('_filtered' if remake else '')
                              + row['#cell_line'] + '@' + row['cells'] + '.tsv')


def open_dfs(df, concat=True, remake=False):
    res_dfs = []
    for index, row in df.iterrows():
        path = get_path(row, remake=remake)
        try:
            tmp_df = pd.read_table(path, comment='#', header=None)
        except EmptyDataError:
            continue
        tmp_df.columns = ['chr', 'pos', 'ref', 'alt', 'BAD'] + ['Q{:.2f}'.format(b) for b in
                                                                [1, 2, 3, 4, 5, 6]] + ['snps', 'cov',
                                                                                       'dataset',
                                                                                       'seg_id',
                                                                                       'p_value']
        group = row['#cell_line'] + '@' + row['cells']
        if not concat:
            tmp_df['cov'] = tmp_df['ref'] + tmp_df['alt']
            tmp_df['max'] = tmp_df[['ref', 'alt']].max(axis=1)
            tmp_df['min'] = tmp_df[['ref', 'alt']].min(axis=1)
            tmp_df['es'] = -np.log2(tmp_df['max'] / tmp_df['min'] / tmp_df['BAD'])
            if res_dfs is None:
                res_dfs = [(group, tmp_df)]
            else:
                res_dfs.append((group, tmp_df))
        else:
            tmp_df['group'] = group
            tmp_df['cell_line'] = row['#cell_line']
            if res_dfs is None:
                res_dfs = [tmp_df]
            else:
                res_dfs.append(tmp_df)
    if concat:
        res_df = pd.concat(res_dfs)
        res_df['cov'] = res_df['ref'] + res_df['alt']
        res_df['max'] = res_df[['ref', 'alt']].max(axis=1)
        res_df['min'] = res_df[['ref', 'alt']].min(axis=1)
        res_df['es'] = -np.log2(res_df['max'] / res_df['min'] / res_df['BAD'])
        return res_df
    else:
        return res_dfs


def update_dist(dist, another_dist):
    for k, v in another_dist.items():
        if k in dist:
            dist[k] += v
        else:
            dist[k] = v
    return dist


def normalize_dist(dist):
    sumr = sum(dist.values())
    return {key: value / sumr for key, value in dist.items()}


def construct_dist_for_cov(cov_df, test_cov_df, metric, norm=False):
    # для данного покрытия строит квантильное распределение
    dist = {}
    for value in cov_df[metric].unique():
        q = percentileofscore(cov_df[metric], value, kind='mean') / 100
        dist[q] = len(test_cov_df[test_cov_df[metric] == value].index)
    if norm:
        dist = normalize_dist(dist)
    return dist


def construct_total_dist(cov_dfs, cov_dfs_test, min_tr, max_tr, metric='es', norm=False):
    total_dist = {}
    covs = []
    total_number_of_snps = sum(len(x.index) for x in cov_dfs_test.values())
    for cov in range(min_tr, max_tr + 1):
        if cov not in cov_dfs:
            continue
        if cov in cov_dfs_test:
            if len(cov_dfs_test[cov].index) >= max(0.0005 * total_number_of_snps, 50):
                covs.append(cov)
                dist = construct_dist_for_cov(cov_dfs[cov],
                                              cov_dfs_test[cov],
                                              metric,
                                              norm=norm)
                total_dist = update_dist(total_dist, dist)
    args = sorted(list(total_dist.keys()))
    vals = [total_dist[x] for x in args]
    return args, vals, covs


def transform_dist_to_list(dist):
    f = lambda *x: sum(x, [])
    return f(*([k] * round(v) for k, v in dist.items()))


def make_dataset(cell_line, lab):
    return '{}@{}'.format(cell_line, lab)


def main(remake=False):
    correlation_file_path = get_correlation_file_path(remake=remake)
    main_df = find_ref_datasets(correlation_file_path)
    res_df = open_dfs(main_df)
    print('DFs concatenated')

    cov_dfs = {cell_line: {} for cell_line in big_cell_lines + ['Other']}
    for cell_line in big_cell_lines + ['Other']:
        if cell_line == 'Other':
            cell_df = res_df[~res_df['cell_line'].isin(big_cell_lines)]
        else:
            cell_df = res_df[res_df['cell_line'] == cell_line]
        for cov in cell_df['cov'].unique():
            cov_dfs[cell_line][cov] = cell_df[cell_df['cov'] == cov].copy()
        print('i split')

    cor_df_test = find_test_datasets(correlation_file_path)
    test_dfs = open_dfs(cor_df_test, concat=False)
    print('Test concatenated')

    min_tr, max_tr = 20, 75

    results = []
    for dataset, dataset_df in test_dfs:
        cov_dfs_test = {}
        for cov in dataset_df['cov'].unique():
            cov_dfs_test[cov] = dataset_df[dataset_df['cov'] == cov].copy()
        print('i split test {}'.format(dataset))
        cell_line = dataset.split('@')[0]
        if cell_line not in big_cell_lines:
            cell_line = 'Other'
        args, vals, covs = construct_total_dist(
            cov_dfs[cell_line],
            cov_dfs_test,
            min_tr=min_tr,
            max_tr=max_tr,
        )
        results.append(
            {
                'args': args,
                'vals': vals,
                'covs': covs,
                'dataset': dataset,
                'snps': len(dataset_df.index)
            })

    cors = pd.read_table(correlation_file_path)

    ref_dists = {x: {} for x in big_cell_lines + ['Other']}
    ref_vars = {x: {} for x in big_cell_lines + ['Other']}

    all_vars = []
    all_metrics = []
    all_cells = []
    all_lines = []
    all_sizes = []

    for d in results:
        if d['args']:
            line, cells = d['dataset'].split('@')
            dist = dict(zip(d['args'], d['vals']))
            cor = cors[
                (cors['#cell_line'] == line)
                & (cors['cells'] == cells)]['cor_by_snp_CAIC'].tolist()
            assert len(cor) == 1
            cor = cor[0]
            if cor > 0.75:
                if line not in big_cell_lines:
                    ref_dists['Other'] = update_dist(ref_dists['Other'], dist)
                else:
                    ref_dists[line] = update_dist(ref_dists[line], dist)

    for key in ref_dists:
        ref_dists[key] = transform_dist_to_list(ref_dists[key])
        ref_vars[key] = np.nanstd(ref_dists[key])

    for d in results:
        if d['args']:
            dist = dict(zip(d['args'], d['vals']))
            line, cells = d['dataset'].split('@')
            snps = d['snps']
            flat_dist = transform_dist_to_list(dist)
            ref_dist = ref_dists[line if line in big_cell_lines else 'Other']
            if not flat_dist or len(flat_dist) == 0:
                continue
            if not ref_dist or len(ref_dist) == 0:
                print('Empty ref dist for {}'.format(line))
                exit(1)
            stat, p = levene(flat_dist, ref_dist)
            all_vars.append((np.nanstd(flat_dist), ref_vars[line if line in big_cell_lines else 'Other']))
            all_metrics.append(p)
            all_cells.append(cells)
            all_lines.append(line)
            all_sizes.append(snps)

    _, all_fdr, _, _ = statsmodels.stats.multitest.multipletests(
        all_metrics, alpha=0.05, method='fdr_bh')

    with open(get_excluded_badmaps_list_path(remake=remake), 'w') as out:
        out.write(pack(['#Cell_line', 'Lab', 'Size', 'dataset_es_var', 'ref_es_var', 'FDR']))
        for fdr, size, line, ce, var in zip(all_fdr, all_sizes, all_lines, all_cells, all_vars):
            out.write(pack([line, ce, size, var[0] ** 2, var[1] ** 2, fdr]))
            # if fdr < 10 ** (-5) and var[0] ** 2 - var[1] ** 2 > 0.05:


if __name__ == '__main__':
    main()
