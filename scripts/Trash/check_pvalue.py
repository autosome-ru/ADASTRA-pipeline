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
    for num in [20]:

        filename = os.path.expanduser('~/fixed_alt_bias_statistics_BAD={:.1f}.tsv'.format(BAD))
        stats = pd.read_table(filename)
        for allele in ('ref', 'alt'):
            stats['{}_counts'.format(allele)] = stats['{}_counts'.format(allele)].astype(int)

        weights_ref = np.load(os.path.expanduser("~/NB_fit_BAD_ref={:.1f}.npy".format(BAD)))
        weights_alt = np.load(os.path.expanduser("~/NB_fit_BAD_alt={:.1f}.npy".format(BAD)))

        scale_df = pd.read_table(os.path.expanduser('~/ref_counts_scaling_BAD={:.1f}.tsv'.format(BAD)))
        scaling = dict(zip(scale_df['allele_reads'], scale_df['new_allele_reads']))

        # alt = []
        # ref = []
        # for ind, row in stats.iterrows():
        #     ref_c, alt_c, counts = row['ref_counts'], row['alt_counts'], row['counts']
        #     print(ref_c, alt_c)
        #     ref_c = scaling[ref_c]
        #     ref += [ref_c] * counts
        #     alt += [alt_c] * counts
        #
        # alt = sorted(alt)
        # ref = sorted(ref)
        #
        # assert len(ref) == len(alt)
        # n = len(ref)
        #
        # plt.scatter(range(n), ref, label='ref')
        # plt.scatter(range(n), alt, label='alt')
        # plt.grid(True)
        # plt.legend()
        # plt.show()

        p_refs = []
        p_alts = []

        max_cover_in_stats = max(
            max(stats['{}_counts'.format('ref' if fixed_allele == 'alt' else 'alt')])
            for fixed_allele in ('ref', 'alt'))

        ref_counts_array = np.zeros(max_cover_in_stats + 1, dtype=np.int64)
        alt_counts_array = np.zeros(max_cover_in_stats + 1, dtype=np.int64)

        for index, row in stats.iterrows():

            ref_c, alt_c, SNP_counts = row['ref_counts'], row['alt_counts'], row['counts']
            ref_c_scaled = round(scaling[ref_c])
            print(ref_c, alt_c)

            if alt_c >= num + d or alt_c < num or weights_ref[alt_c][0] == 0:
                p_ref = 1
            else:
                p_ref = calculate_p_value(ref_c_scaled, weights_ref[alt_c][1], weights_ref[alt_c][0], 1 / (BAD + 1))

            if ref_c >= num + d or ref_c < num or weights_alt[ref_c][0] == 0:
                p_alt = 1
            else:
                p_alt = calculate_p_value(alt_c, weights_alt[ref_c][1], weights_alt[ref_c][0], 1 / (BAD + 1))

            print(p_ref, p_alt)
            p_refs.append(p_ref)
            p_alts.append(p_alt)

            if p_ref <= 1/(10**5):
                ref_raw = scaling[ref_c]
                ref_floor = np.floor(ref_raw)
                ref_ceil = np.ceil(ref_raw)
                part = ref_raw - ref_floor
                floor_counts = np.ceil(SNP_counts * (1 - part))
                ceil_counts = SNP_counts - floor_counts

                # print(k_raw, k_floor, k_ceil, SNP_counts, floor_counts, ceil_counts, part)

                # counts_array[int(k_floor)] += SNP_counts
                ref_counts_array[int(ref_floor)] += floor_counts
                ref_counts_array[int(ref_ceil)] += ceil_counts
            elif p_alt <= 1/(10**5):

                alt_counts_array[alt_c] += SNP_counts


        stats['ref_p'] = p_refs
        stats['alt_p'] = p_alts

        fig, ax = plt.subplots(figsize=(10, 8))
        x = np.array(range(1, 18))

        palts = stats[stats['alt_p'] != 1.0]
        sns.barplot(x=x, y=[palts[(palts['alt_p'] <= 1 / 10 ** k) & (palts['alt_p'] > 1 / 10 ** (k + 1))]['counts'].sum() for k in x], label='alt', ax=ax,
                    color='C1', alpha=0.5)

        prefs = stats[stats['ref_p'] != 1.0]
        sns.barplot(x=x, y=[prefs[(prefs['ref_p'] <= 1 / 10 ** k) & (prefs['ref_p'] > 1 / 10 ** (k + 1))]['counts'].sum() for k in x], label='ref', ax=ax,
                    color='C0', alpha=0.5)

        plt.grid(True)
        plt.legend()
        plt.title('ref-alt p_value on BAD={:.1f}\n counts {} - {}'.format(BAD, num, num+d))
        plt.xlabel('x: -log10 p_value = x')
        plt.ylabel('snp count')
        plt.savefig(os.path.expanduser('~/fixed_alt/p_dist_BAD={:.1f}_{}_{}.png'.format(BAD, num, num+d)))
        plt.show()
        plt.close(fig)

        ref_df = pd.DataFrame({'counts': ref_counts_array, 'type': 'ref'})

        alt_df = pd.DataFrame({'counts': alt_counts_array, 'type': 'alt'})

        plot_df = ref_df.append(alt_df)

        fig, ax = plt.subplots(figsize=(10, 8))

        try:

            plt.scatter(x=alt_df.index,
                        y=(alt_df['counts'])/alt_df['counts'].sum(), label='alt',
                        color='C1', alpha=0.5)

            plt.scatter(x=ref_df.index,
                        y=(ref_df['counts'])/ref_df['counts'].sum(), label='ref',
                        color='C0', alpha=0.5)
        except KeyboardInterrupt:
            print("I'm busy plotting")

        plt.grid(True)
        plt.legend()
        plt.title('density {} - {}'.format(num, num+d))
        plt.xlabel('read_counts')
        plt.ylabel('snp counts')
        plt.savefig(os.path.expanduser('~/fixed_alt/net_distributions_BAD={:.1f}_{}_{}.png'.format(BAD, num, num+d)))
        plt.show()
