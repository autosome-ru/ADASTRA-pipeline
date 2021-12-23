import pandas as pd
import os
import numpy as np

from scripts.HELPERS.helpers import get_states
from scripts.HELPERS.paths import get_heatmap_data_path

cell_sign_dict = {
    'K562': 'K562__myelogenous_leukemia_',
    'MCF7': "MCF7__Invasive_ductal_breast_carcinoma_",
    'A549': 'A549__lung_carcinoma_',
    'HCT116': 'HCT-116__colon_carcinoma_',
    '22RV1': '22RV1__prostate_carcinoma_',
}


def cell_line_in_file_from_sign(cell_sign, file):
    check_other = True
    for sign, cell_line in cell_sign_dict.items():
        if cell_line in file:
            if cell_sign == sign:
                return True
            check_other = False
    if check_other and cell_sign == 'Other':
        return True
    else:
        return False


def main(states_sign, b_penalty):
    model = 'CAIC@{}@{}'.format(states_sign, b_penalty)
    print(model)
    for cell_sign in ('K562', 'MCF7', 'A549', 'HCT116', '22RV1', 'Other'):
        print(cell_sign)
        states = get_states(states_sign)
        heatmap_dir = os.path.join(get_heatmap_data_path(), model + '_tables/')
        dfs = []
        for file in os.listdir(heatmap_dir):
            if not cell_line_in_file_from_sign(cell_sign, file):
                continue
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

            res_dir = os.path.expanduser('~/PARAMETERS/counts/')
            sum_df.to_csv(os.path.join(res_dir, 'counts_deltaqm_{}_{}_{:.2f}.tsv'.format(cell_sign, model, BAD)),
                          index=False, sep='\t')