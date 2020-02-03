import pandas as pd
import os
import sys
import numpy as np

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.helpers import states

for BAD in states:
    df = pd.read_table(os.path.expanduser('~/unionSNPs.tsv'))
    df.columns = ['chr', 'pos', 'cov', 'BAD', 'COSMIC'] + ['Q{:.2f}'.format(state) for state in states]
    df['deltaQ{:.2f}'.format(BAD)] = df['Q{:.2f}'.format(BAD)] - df[
        ['Q{:.2f}'.format(another_BAD) for another_BAD in states if another_BAD != BAD]].max(axis=1)
    print(df['deltaQ{:.2f}'.format(BAD)].unique())
    min_tr = df['deltaQ{:.2f}'.format(BAD)].min()
    max_tr = df['deltaQ{:.2f}'.format(BAD)].max()
    thresholds = list(np.linspace(min_tr, -100, min(50, min_tr + 100))) if min_tr < -100 else [] + \
        list(range(-100, 101, 1)) + list(np.linspace(101, max_tr, min(50, max_tr - 100))) if max_tr > 100 else []
    print(thresholds)
    sum_df = None
    for threshold in thresholds:
        print('Now doing threshold = {}, BAD={:.2f}'.format(threshold, BAD))
        df = df[df['deltaQ{:.2f}'.format(BAD)] >= threshold]
        df_counts = df.groupby(['BAD', 'COSMIC']).size().reset_index(name='counts')
        df_counts = df_counts.groupby(['BAD', 'COSMIC'], as_index=False)['counts'].sum()
        df_counts['threshold'] = threshold
        if sum_df is None:
            sum_df = df_counts
        else:
            sum_df = sum_df.append(df_counts)
    sum_df.to_csv(os.path.expanduser('~/PARAMETERS/counts/counts_deltaqm_{:.2f}.tsv'.format(BAD)), index=False, sep='\t')
