import os
import numpy as np
import pandas as pd

from scripts.HELPERS.helpers import states
from scripts.HELPERS.paths import create_neg_bin_weights_path_function, create_neg_bin_stats_path_function, \
    get_release_stats_path


def main():
    column_names = ['r', 'w', 'status', 'gof']
    for BAD in states:
        counts_df = pd.read_table(create_neg_bin_stats_path_function(BAD))
        for allele in ('ref', 'alt'):
            np_weights = np.load(create_neg_bin_weights_path_function(allele, BAD))
            df = pd.DataFrame(columns=column_names)
            for i in range(len(column_names)):
                df[column_names[i]] = np_weights[:, i]
            counts = []
            for fix_c in df.index:
                counts.append(counts_df[counts_df['{}_counts'.format(allele)] == fix_c]['counts'].sum())
            df['allele_reads'] = counts
            df.to_csv(os.path.join(get_release_stats_path(), 'NBweights_{}_BAD{:.1f}.tsv'.format(allele, BAD)),
                      index=False, sep='\t')
