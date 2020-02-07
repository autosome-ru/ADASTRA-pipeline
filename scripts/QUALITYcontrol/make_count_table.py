import pandas as pd
import os
import sys
import numpy as np

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.helpers import states

for BAD in states:
    df = pd.read_table(os.path.expanduser('~/unionSNPs.tsv'))
    df.columns = ['chr', 'pos', 'cov', 'BAD', 'COSMIC'] + ['Q{:.2f}'.format(state) for state in states]
    df['BAD'] = BAD
    df['threshold'] = df['Q{:.2f}'.format(BAD)] - df[
        ['Q{:.2f}'.format(another_BAD) for another_BAD in states if another_BAD != BAD]].max(axis=1)
    print(df['threshold'].unique())
    min_tr = df['threshold'].min()
    max_tr = df['threshold'].max()
    N = 300
    idxs = set(int(x) for x in np.linspace(0, len(df.index) - 1, N))
    thresholds = []
    for index, row in df.sort_values(by='threshold').reset_index(drop=True).iterrows():
        if index in idxs:
            thresholds.append(row['threshold'])
    print(thresholds)
    sum_df = None
    for threshold in thresholds:
        print('Now doing threshold = {}, BAD={:.2f}'.format(threshold, BAD))
        print('before: {}'.format(len(df.index)))
        df = df[df['threshold'] >= threshold]
        print('after: {}'.format(len(df.index)))
        df_counts = df.groupby(['BAD', 'COSMIC']).size().reset_index(name='counts')
        df_counts = df_counts.groupby(['BAD', 'COSMIC'], as_index=False)['counts'].sum()
        df_counts['threshold'] = threshold
        if sum_df is None:
            sum_df = df_counts
        else:
            sum_df = sum_df.append(df_counts)
    sum_df.to_csv(os.path.expanduser('~/PARAMETERS/counts/counts_deltaqm_{:.2f}.tsv'.format(BAD)), index=False,
                  sep='\t')
