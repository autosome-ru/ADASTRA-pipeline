import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.FITnoise.true_quantile_norm import get_normalized_allele_reads

for flag in ('q_norm', 'change'):
    for BAD in [1, 2, 3, 4, 5, 6, 4 / 3, 5 / 2, 3 / 2]:
        stats = pd.read_table(os.path.expanduser('~/bias_statistics_BAD={:.1f}.tsv'.format(BAD)))
        #stats = stats[stats['allele_reads'] <= 2000]
        fig, ax = plt.subplots(figsize=(10, 8))

        assert list(stats['allele_reads']) == sorted(list(stats['allele_reads']))
        stats['new_allele_reads'] = get_normalized_allele_reads(stats)
        #stats = stats.fillna(0)
        #print(stats[['allele_reads', 'new_allele_reads']])

        if flag == 'q_norm':
            plt.scatter(x=stats['allele_reads'], y=stats['new_allele_reads'], color='C0', label='ref_norm', s=50)
            plt.plot([0, max(stats['allele_reads'])], [0, max(stats['allele_reads'])], label='y=x', color='black')

            plt.xlabel('old_ref')
            plt.ylabel('new_ref')
            plt.title('ref read counts before and after rank normalization, BAD={}'.format(BAD))

        if flag == 'change':
            sns.scatterplot(x=stats['allele_reads'], y=np.log10(stats['ref'] + 1), color='C3', ax=ax, label='ref',
                            s=50)
            sns.scatterplot(x=stats['allele_reads'], y=np.log10(stats['alt'] + 1), color='C1', ax=ax, label='alt', s=50)
            sns.scatterplot(x=stats['new_allele_reads'], y=np.log10(stats['ref'] + 1), color='C0', ax=ax, label='ref_norm',
                            s=50)
            plt.xlabel('allele read counts')
            plt.ylabel('log10 snp counts')
            plt.title('before and after REF rank normalization, BAD={}'.format(BAD))
        if flag in ('ref', 'before'):
            sns.scatterplot(x=stats['allele_reads'], y=stats['ref'], color='C0', ax=ax, label='ref', s=50)
        if flag in ('alt', 'before'):
            sns.scatterplot(x=stats['allele_reads'], y=stats['alt'], color='C1', ax=ax, label='alt', s=40)
        if flag == 'ratio':
            sns.scatterplot(x=stats['allele_reads'], y=np.log10(stats['ref'] + 1) - np.log10(1 + stats['alt']), ax=ax,
                            color='C5', label='log10 (ref+1)/(alt+1)', s=40)
            plt.axhline(y=0, color='black')

        plt.grid(True)
        # if flag == 'ratio':
        #     plt.title('ref-alt ratio for BAD={:.1f}'.format(BAD))
        #     plt.ylabel('log10 (ref+1)/(alt+1)')
        # else:
        #     plt.title('ref-alt bias for BAD={:.1f}'.format(BAD))
            #plt.ylabel('log10 count')
        plt.legend()
        plt.savefig(os.path.expanduser('~/qnorm/bias_q_{}_BAD={:.1f}.png'.format(flag, BAD)))
        plt.close(fig)
