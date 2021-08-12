import json

import pandas as pd
from pandas.errors import EmptyDataError
import os
import numpy as np
from scipy.stats import levene
import statsmodels.stats.multitest

from scripts.HELPERS.helpers import pack, segmentation_states
from scripts.HELPERS.paths import get_excluded_badmaps_list_path, get_correlation_file_path, get_correlation_path


def find_test_datasets(cor_stats_path):
    df = pd.read_table(cor_stats_path)
    return df


def get_path(row, remake=False):
    return os.path.join(get_correlation_path(), 'CAIC_tables{}/'.format('_filtered' if remake else '')
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
                                                               segmentation_states] + ['snps', 'cov',
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


def construct_dist_for_cov(test_cov_df):
    dist = {}
    for value in test_cov_df['es'].unique():
        dist[value] = len(test_cov_df[test_cov_df['es'] == value].index)
    return dist


def construct_total_dist(cov_dfs_test, min_tr, max_tr):
    total_dist = {}
    covs = []
    for cov in range(min_tr, max_tr + 1):
        if cov in cov_dfs_test:
            if len(cov_dfs_test[cov].index) >= 50:
                covs.append(cov)
                dist = construct_dist_for_cov(cov_dfs_test[cov])
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
    cor_df_test = find_test_datasets(correlation_file_path)
    test_dfs = open_dfs(cor_df_test, remake=remake, concat=False)
    print('Test concatenated')

    min_tr, max_tr = 20, 75

    results = []
    for dataset, dataset_df in test_dfs:
        cov_dfs_test = {}
        for cov in dataset_df['cov'].unique():
            cov_dfs_test[cov] = dataset_df[dataset_df['cov'] == cov].copy()
        print('Split test {}'.format(dataset))
        args, vals, covs = construct_total_dist(cov_dfs_test, min_tr=min_tr, max_tr=max_tr)
        results.append({'args': args, 'vals': vals, 'covs': covs, 'dataset': dataset, 'snps': len(dataset_df.index)})

    with open(os.path.expanduser('~/cov_res_debug.json'), 'w') as f:
        json.dump(results, f)

    cors = pd.read_table(correlation_file_path)

    if not remake:
        # collect_stats
        cell_line_data = {}
        for d in results:
            if d['args']:
                line, cells = d['dataset'].split('@')
                cor = cors[(cors['#cell_line'] == line) & (cors['cells'] == cells)]['cor_by_snp_CAIC'].tolist()
                assert len(cor) == 1
                cor = cor[0]
                if not pd.isna(cor):
                    cell_line_data.setdefault(line, {
                        'correlations': [],
                        'cells': [],
                        'snps': []
                    })
                    cell_line_data[line]['correlations'].append(cor)
                    cell_line_data[line]['cells'].append(cells)
                    cell_line_data[line]['snps'].append(sum(d['vals']))

        # construct big cell lines
        big_cell_lines = set()
        cell_line_reference = {}
        for line, data in cell_line_data.items():
            if len(data['correlations']) < 4:
                continue
            cor_treshold = np.quantile(data['correlations'], 0.75)
            datasets = [(cells, cor, snps) for cells, cor, snps in
                        zip(data['cells'], data['correlations'], data['snps']) if cor >= cor_treshold]
            snps = sum(x[2] for x in datasets)
            if snps >= 25000:
                cell_line_reference[line] = [x[0] for x in datasets]
                big_cell_lines.add(line)
    else:
        prev_excluded = pd.read_table(get_excluded_badmaps_list_path(remake=False))
        big_cell_lines = set()
        cell_line_reference = {}
        for index, row in prev_excluded.iterrows():
            if row['is_ref']:
                big_cell_lines.add(row['#cell_line'])
                cell_line_reference.setdefault(row['#cell_line'], []).append(row['sample'])

    big_cell_lines = list(big_cell_lines)

    ref_dists = {x: {} for x in big_cell_lines + ['Other']}
    ref_vars = {x: {} for x in big_cell_lines + ['Other']}

    all_vars = []
    all_metrics = []
    all_cells = []
    all_lines = []
    all_sizes = []
    all_is_ref = []

    for d in results:
        if d['args']:
            line, cells = d['dataset'].split('@')
            dist = dict(zip(d['args'], d['vals']))
            if line in big_cell_lines:
                if cells in cell_line_reference[line]:
                    ref_dists[line] = update_dist(ref_dists[line], dist)
                    ref_dists['Other'] = update_dist(ref_dists[line], dist)

    for key in ref_dists:
        ref_dists[key] = transform_dist_to_list(ref_dists[key])
        ref_vars[key] = np.nanstd(ref_dists[key])

    print(ref_vars)

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
            assert not pd.isna(p)
            all_vars.append((np.nanstd(flat_dist), ref_vars[line if line in big_cell_lines else 'Other']))
            all_metrics.append(p)
            all_cells.append(cells)
            all_lines.append(line)
            all_sizes.append(snps)
            all_is_ref.append(True if line in big_cell_lines and cells in cell_line_reference[line] else False)

    _, all_fdr, _, _ = statsmodels.stats.multitest.multipletests(
        all_metrics, alpha=0.05, method='fdr_bh')

    with open(get_excluded_badmaps_list_path(remake=remake), 'w') as out:
        out.write(pack(['#cell_line', 'sample', 'size', 'dataset_es_var', 'ref_es_var', 'fdr', 'is_ref']))
        for fdr, size, line, ce, var, ref in zip(all_fdr, all_sizes, all_lines, all_cells, all_vars, all_is_ref):
            out.write(pack([line, ce, size, var[0] ** 2, var[1] ** 2, fdr, ref]))
            # if fdr < 10 ** (-50) and var[0] ** 2 - var[1] ** 2 > 0.1 and not ref:


if __name__ == '__main__':
    main()
