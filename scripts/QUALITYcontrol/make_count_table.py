import pandas as pd
import os

df = pd.read_table(os.path.expanduser('~/uniuonSNPs.tsv'))
df.columns = ['chr', 'pos', 'cov', 'BAD', 'COSMIC']
thresholds = [x for x in range(10, 51, 5)]
for threshold in thresholds:
    print('Now doing threshold = {}'.format(threshold))
    df = df[df['cov'] > threshold]
    df_counts = df.groupby(['BAD', 'COSMIC']).size().reset_index(name='counts')
    df_counts = df_counts.groupby(['BAD', 'COSMIC'], as_index=False)['counts'].sum()
    df_counts.to_csv(os.path.expanduser(
        '~/PARAMETERS/counts/counts_covm_{}.tsv'.format(threshold)), index=False, sep='\t')
