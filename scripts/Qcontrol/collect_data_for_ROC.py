import pandas as pd
import os
import numpy as np

from scripts.HELPERS.helpers import get_states, get_params_from_model_name
from scripts.HELPERS.paths import get_heatmap_data_path, get_release_stats_path


def main(model):
    print(model)
    states_sign = get_params_from_model_name(model)['states_set']
    states = get_states(states_sign)
    heatmap_dir = os.path.join(get_heatmap_data_path(), model + '_tables/')
    dfs = []
    for file in os.listdir(heatmap_dir):
        try:
            dfs.append(pd.read_table(os.path.join(heatmap_dir, file), header=None, comment='#'))
        except pd.errors.EmptyDataError:
            continue
    full_df = pd.concat(dfs)
    full_df.columns = ['chr', 'pos', 'cov', 'BAD', 'COSMIC'] + ['Q{:.2f}'.format(state) for state in states]
    print('i read')

    for BAD in states:
        df = full_df.copy()
        df['BAD'] = BAD
        df['threshold'] = df['Q{:.2f}'.format(BAD)] - df[
            ['Q{:.2f}'.format(another_BAD) for another_BAD in states if another_BAD != BAD]].max(axis=1)
        print(df['threshold'].unique())
        df = df[['BAD', 'COSMIC', 'threshold']]
        # min_tr = df['threshold'].min()
        # max_tr = df['threshold'].max()
        N = 1000
        idxs = set(int(x) for x in np.linspace(0, len(df.index) - 1, N))
        sorted_by_thresholds = df['threshold'].to_numpy(copy=True)
        sorted_by_thresholds.sort()
        print('isort')
        thresholds = {0}
        for index, value in enumerate(sorted_by_thresholds):
            if index in idxs:
                thresholds.add(value)
        thresholds = sorted(list(thresholds))
        print(thresholds)
        sum_df = None
        for threshold in thresholds:
            print('Now doing threshold = {}, BAD={:.2f}'.format(threshold, BAD))
            before = len(df.index)
            df = df[df['threshold'] >= threshold]
            print('change: {}'.format(before - len(df.index)))
            df_counts = df.groupby(['BAD', 'COSMIC']).size().reset_index(name='counts')
            df_counts = df_counts.groupby(['BAD', 'COSMIC'], as_index=False)['counts'].sum()
            df_counts['threshold'] = threshold
            if sum_df is None:
                sum_df = df_counts
            else:
                sum_df = sum_df.append(df_counts)

        res_dir = os.path.join(get_release_stats_path(), 'roc_data')
        sum_df.to_csv(os.path.join(res_dir, 'counts_deltaqm_{}_{:.2f}.tsv'.format(model, BAD)), index=False, sep='\t')
