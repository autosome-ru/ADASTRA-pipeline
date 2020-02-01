import pandas as pd
import os

df = pd.read_table(os.path.expanduser('~/uniuonSNPs.tsv'))
df.columns = ['chr', 'pos', 'cov', 'BAD', 'COSMIC']
thresholds = [x for x in range(10, 101, 1)]
sum_df = None
for threshold in thresholds:
    print('Now doing threshold = {}'.format(threshold))
    df = df[df['cov'] > threshold]
    df_counts = df.groupby(['BAD', 'COSMIC']).size().reset_index(name='counts')
    df_counts = df_counts.groupby(['BAD', 'COSMIC'], as_index=False)['counts'].sum()
    df_counts['threshold'] = threshold
    if sum_df is None:
        sum_df = df_counts
    else:
        sum_df.append(df_counts)
sum_df.to_csv(os.path.expanduser('~/PARAMETERS/counts/counts_m.tsv'), index=False, sep='\t')
