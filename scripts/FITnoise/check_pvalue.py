import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from scipy import stats as st


def calculate_p_value(main_c, w, r, p):
    assert 0 < w <= 1
    assert r > 1
    print(r, w)
    dist1 = st.nbinom(r, p)
    cdf1 = dist1.cdf
    dist2 = st.nbinom(r, 1 - p)
    cdf2 = dist2.cdf
    return (1 - cdf1(main_c)) * w + \
           (1 - cdf2(main_c)) * (1 - w)


for BAD in [1]:

    d = 10
    for num in [5]:

        filename = os.path.expanduser('~/fixed_alt_bias_statistics_BAD={:.1f}.tsv'.format(BAD))
        stats = pd.read_table(filename)
        for allele in ('ref', 'alt'):
            stats['{}_counts'.format(allele)] = stats['{}_counts'.format(allele)].astype(int)

        weights = np.load(os.path.expanduser("~/NB_fit_BAD={:.1f}.npy".format(BAD)))

        scale_df = pd.read_table(os.path.expanduser('~/ref_counts_scaling_BAD={:.1f}.tsv'.format(BAD)))
        scaling = dict(zip(scale_df['allele_reads'], scale_df['new_allele_reads']))

        p_refs = []
        p_alts = []
        for index, row in stats.iterrows():

            ref_c, alt_c = row['ref_counts'], row['alt_counts']
            ref_c = int(scaling[ref_c])
            print(ref_c, alt_c)

            if alt_c >= num + d or alt_c < num or weights[alt_c][0] == 0:
                p_ref = 1
            else:
                p_ref = calculate_p_value(ref_c, weights[alt_c][1], weights[alt_c][0], 1 / (BAD + 1))

            if ref_c >= num + d or alt_c < num or weights[ref_c][0] == 0:
                p_alt = 1
            else:
                p_alt = calculate_p_value(alt_c, weights[ref_c][1], weights[ref_c][0], 1 / (BAD + 1))

            print(p_ref, p_alt)
            p_refs.append(p_ref)
            p_alts.append(p_alt)

        stats['ref_p'] = p_refs
        stats['alt_p'] = p_alts

        fig, ax = plt.subplots(figsize=(10, 8))
        x = np.array(range(5, 18))

        palts = stats[stats['alt_p'] != 1.0]
        sns.barplot(x=x, y=[palts[palts['alt_p'] <= 1 / 10 ** k]['counts'].sum() for k in x], label='alt', ax=ax,
                    color='C1', alpha=0.5)

        prefs = stats[stats['ref_p'] != 1.0]
        sns.barplot(x=x, y=[prefs[prefs['ref_p'] <= 1 / 10 ** k]['counts'].sum() for k in x], label='ref', ax=ax,
                    color='C0', alpha=0.5)

        plt.grid(True)
        plt.legend()
        plt.title('ref-alt p_value on BAD={:.1f}\n counts {} - {}'.format(BAD, num, num+d))
        plt.xlabel('x: -log10 p_value >= x')
        plt.ylabel('snp count')
        plt.savefig(os.path.expanduser('~/fixed_alt/p_dist_BAD={:.1f}_{}_{}.png'.format(BAD, num, num+d)))
        plt.show()
