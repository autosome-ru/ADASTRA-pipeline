import pandas as pd
import os
import sys
import numpy as np

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.helpers import states

for BAD in states:
    df = pd.read_table(os.path.expanduser('~/unionSNPs_{:.2f}.tsv'.format(BAD)))
    df.columns = ['chr', 'pos', 'cov', 'BAD', 'COSMIC'] + ['Q{:.2f}'.format(state) for state in states]
    df['BAD'] = BAD
    df['threshold'] = df['Q{:.2f}'.format(BAD)]
    print(df['threshold'].unique())
    min_tr = df['threshold'].min()
    max_tr = df['threshold'].max()
    N = 300
    idxs = set(int(x) for x in np.linspace(0, len(df.index) - 1, N))
    thresholds = []
    sorted_by_thresholds = df['threshold'].to_numpy(copy=True)
    sorted_by_thresholds.sort()
    print('isort')
    for index, value in enumerate(sorted_by_thresholds):
        if index in idxs:
            thresholds.append(value)
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
    sum_df.to_csv(os.path.expanduser('~/PARAMETERS/counts/counts_rawlm_{:.2f}.tsv'.format(BAD)), index=False,
                  sep='\t')
