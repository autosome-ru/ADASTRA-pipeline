import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from scipy import stats as st

filename = os.path.expanduser('~/CTCF_HUMAN_SNP_table.tsv')
sum_df = pd.read_table(filename)
BAD = 1
sum_df = sum_df[sum_df['BAD'] == BAD]
alt_counts = sorted(sum_df['alt_read_counts'].tolist())
sum_df = sum_df.sort_values('ref_read_counts')
sum_df.reset_index()
sum_df['ref_read_counts'] = alt_counts

out_t = pd.DataFrame()
out_t['alt_counts'] = sum_df['alt_read_counts']
out_t['ref_counts'] = sum_df['ref_read_counts']
out_t = out_t.groupby(['alt_counts', 'ref_counts']).size().reset_index(name='counts')
out_t.fillna(0, inplace=True)
print(out_t)
with open(os.path.expanduser('~/fixed_alt_bias_statistics_BAD={:.1f}.tsv'.format(BAD)), 'w') as out:
    out_t.to_csv(out, sep="\t", index=False)
