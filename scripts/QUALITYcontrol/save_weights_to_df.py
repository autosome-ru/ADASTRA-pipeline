import sys
import numpy as np
import pandas as pd

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths_for_components import parameters_path
from scripts.HELPERS.helpers import states


column_names = ['r', 'w', 'status', 'gof']
for BAD in states:
    ref = np.load(parameters_path + 'NBweights_ref_BAD={:.1f}.npy'.format(BAD))
    alt = np.load(parameters_path + 'NBweights_alt_BAD={:.1f}.npy'.format(BAD))
    counts_df = pd.read_table(parameters_path + 'fixed_alt_bias_statistics_BAD={:.1f}.tsv'.format(BAD))
    ref_df = pd.DataFrame(columns=column_names)
    alt_df = pd.DataFrame(columns=column_names)
    for i in range(len(column_names)):
        ref_df[column_names[i]] = ref[:, i]
        alt_df[column_names[i]] = alt[:, i]
    ref_counts = []
    alt_counts = []
    for fix_c in ref_df.index:
        ref_counts.append(counts_df[counts_df['ref_counts'] == fix_c]['counts'].sum())
        alt_counts.append(counts_df[counts_df['alt_counts'] == fix_c]['counts'].sum())
    ref_df['allele_reads'] = ref_counts
    alt_df['allele_reads'] = alt_counts
    ref_df.to_csv(parameters_path + 'weights/NBweights_ref_BAD={:.1f}.tsv'.format(BAD), sep='\t')
    alt_df.to_csv(parameters_path + 'weights/NBweights_alt_BAD={:.1f}.tsv'.format(BAD), sep='\t')
