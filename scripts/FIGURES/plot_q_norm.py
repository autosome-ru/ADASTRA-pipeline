import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.FITnoise.quantile_norm import normalize_values

for flag in ('change',):
    for BAD in [1, 2, 3, 4, 5, 6, 4 / 3, 5 / 2, 3 / 2]:
        stats = pd.read_table(os.path.expanduser('~/bias_statistics_BAD={:.1f}.tsv'.format(BAD)))
        #stats = stats[stats['allele_reads'] <= 500]
        fig, ax = plt.subplots(figsize=(10, 8))

        values = stats.index
        ref_counts = stats['ref']
        alt_counts = stats['alt']
        stats['new_allele_reads'] = normalize_values(values, values, ref_counts, alt_counts)

        if flag == 'q_norm':
            sns.scatterplot(x=stats['new_allele_reads'], y=np.log10(stats['ref'] + 1), color='C0', ax=ax, label='ref_norm', s=50)
            sns.scatterplot(x=stats['allele_reads'], y=np.log10(stats['alt'] + 1), color='C1', ax=ax, label='alt', s=50)
        if flag == 'change':
            sns.scatterplot(x=stats['allele_reads'], y=np.log10(stats['ref'] + 1), color='C3', ax=ax, label='ref',
                            s=50)
            sns.scatterplot(x=stats['allele_reads'], y=np.log10(stats['alt'] + 1), color='C1', ax=ax, label='alt', s=50)
            sns.scatterplot(x=stats['new_allele_reads'], y=np.log10(stats['ref'] + 1), color='C0', ax=ax, label='ref_norm',
                            s=50)
        if flag in ('ref', 'before'):
            sns.scatterplot(x=stats['allele_reads'], y=stats['ref'], color='C0', ax=ax, label='ref', s=50)
        if flag in ('alt', 'before'):
            sns.scatterplot(x=stats['allele_reads'], y=stats['alt'], color='C1', ax=ax, label='alt', s=40)
        if flag == 'ratio':
            sns.scatterplot(x=stats['allele_reads'], y=np.log10(stats['ref'] + 1) - np.log10(1 + stats['alt']), ax=ax,
                            color='C5', label='log10 (ref+1)/(alt+1)', s=40)
            plt.axhline(y=0, color='black')

        plt.grid(True)
        if flag == 'ratio':
            plt.title('ref-alt ratio for BAD={:.1f}'.format(BAD))
            plt.ylabel('log10 (ref+1)/(alt+1)')
        else:
            plt.title('ref-alt bias for BAD={:.1f}'.format(BAD))
            plt.ylabel('log10 count')
        plt.legend()
        plt.savefig(os.path.expanduser('~/qnorm/bias_q_{}_BAD={:.1f}.png'.format(flag, BAD)))
        plt.close(fig)
