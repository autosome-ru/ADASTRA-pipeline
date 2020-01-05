import os
import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, ticker
from scipy import optimize
from scipy import stats as st
from scipy.special import beta, psi, poch, comb
from sklearn import metrics
import seaborn as sns
import subprocess

for BAD in [1]:
    filename = os.path.expanduser('~/fixed_alt_bias_statistics_BAD={:.1f}.tsv'.format(BAD))
    stats = pd.read_table(filename)
    for allele in ('ref', 'alt'):
        stats['{}_counts'.format(allele)] = stats['{}_counts'.format(allele)].astype(int)

    scale_df = pd.read_table(os.path.expanduser('~/ref_counts_scaling_BAD={:.1f}.tsv'.format(BAD)))
    scaling = dict(zip(scale_df['allele_reads'], scale_df['new_allele_reads']))

    max_cover_in_stats = max(
        max(stats['{}_counts'.format('ref' if fixed_allele == 'alt' else 'alt')])
        for fixed_allele in ('ref', 'alt'))

    ref_counts_array = np.zeros(max_cover_in_stats + 1, dtype=np.int64)
    alt_counts_array = np.zeros(max_cover_in_stats + 1, dtype=np.int64)

    print(max_cover_in_stats)

    for index, row in stats.iterrows():
        ref, alt, SNP_counts = row['ref_counts'], row['alt_counts'], row['counts']
        print(ref, alt)
        ref_raw = scaling[ref]
        ref_floor = np.floor(ref_raw)
        ref_ceil = np.ceil(ref_raw)
        part = ref_raw - ref_floor
        floor_counts = np.ceil(SNP_counts * (1 - part))
        ceil_counts = SNP_counts - floor_counts

        # print(k_raw, k_floor, k_ceil, SNP_counts, floor_counts, ceil_counts, part)

        # counts_array[int(k_floor)] += SNP_counts
        ref_counts_array[int(ref_floor)] += floor_counts
        ref_counts_array[int(ref_ceil)] += ceil_counts

        alt_counts_array[alt] += SNP_counts

    ref_df = pd.DataFrame({'counts': ref_counts_array, 'type': 'ref'})

    alt_df = pd.DataFrame({'counts': alt_counts_array, 'type': 'alt'})

    plot_df = ref_df.append(alt_df)

    fig, ax = plt.subplots(figsize=(10, 8))

    try:

        plt.scatter(x=alt_df.index,
                    y=np.log10(alt_df['counts']), label='alt',
                    color='C1', alpha=0.5)

        plt.scatter(x=ref_df.index,
                    y=np.log10(ref_df['counts']), label='ref',
                    color='C0', alpha=0.5)
    except KeyboardInterrupt:
        print("I'm busy plotting")


    plt.grid(True)
    plt.legend()
    plt.title('Log density')
    plt.xlabel('read_counts')
    plt.ylabel('log10 snp counts')
    plt.savefig(os.path.expanduser('~/fixed_alt/net_distributions_BAD={}.png'.format(BAD)))
    plt.show()
