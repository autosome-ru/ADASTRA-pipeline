import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, ticker
import seaborn as sns
from scipy import stats as st

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.helpers import states, read_weights


def make_negative_binom_density(r, p, w, size_of_counts):
    negative_binom_density_array = np.zeros(size_of_counts + 1, dtype=np.float128)
    dist1 = st.nbinom(r, p)
    f1 = dist1.pmf
    cdf1 = dist1.cdf
    dist2 = st.nbinom(r, 1 - p)
    f2 = dist2.pmf
    cdf2 = dist2.cdf
    negative_binom_norm = (cdf1(size_of_counts) - cdf1(4)) * w + \
                          (cdf2(size_of_counts) - cdf2(4)) * (1 - w)

    for k in range(5, size_of_counts + 1):
        negative_binom_density_array[k] = (w * f1(k) + (1 - w) * f2(k)) / negative_binom_norm
    return negative_binom_density_array


r_dict, w_dict, gof_dict = read_weights()
fix_c = 20
main_allele = "alt"
fixed_allele = "ref" if main_allele == "alt" else "alt"
for BAD in states:
    stats = pd.read_table()
    stats_filtered = stats[stats['{}_counts'.format(fixed_allele)] == fix_c]
    max_cover_in_stats = max(stats_filtered['{}_counts'.format(main_allele)])
    counts_array = np.zeros(max_cover_in_stats + 1, dtype=np.int64)
    for index, row in stats_filtered.iterrows():
        k, SNP_counts = row['{}_counts'.format(main_allele)], row['counts']
        counts_array[k] = SNP_counts

    total_snps = counts_array[0:max_cover_in_stats + 1].sum()
    fig, ax = plt.subplots(figsize=(10, 8))
    x = list(range(max_cover_in_stats + 1))
    sns.barplot(x=x,
                y=counts_array[0:max_cover_in_stats + 1] / total_snps, ax=ax)
    ax.xaxis.set_major_locator(ticker.FixedLocator(np.arange(0, len(x), 5)))
    ax.xaxis.set_major_formatter(ticker.FixedFormatter(x[::5]))
    ax.tick_params(axis="x", rotation=90)

    r, w, gof = (r_dict[fixed_allele][BAD][fix_c],
                 w_dict[fixed_allele][BAD][fix_c],
                 gof_dict[fixed_allele][BAD][fix_c])
    current_density = make_negative_binom_density(r, BAD, w, fix_c)

    asb = counts_array[5:].sum()

    label = 'negative binom fit for {}r_dict, w_dict, ' \
            'gof_dict = read_weights()\ntotal observations: {}\nr={:.2f}, p={:.2f}, q={}, w={:.2f}\ngof={:.4}'.format(
        main_allele, total_snps, r, BAD, w, asb, gof)
    plt.plot(list(range(fix_c + 1)), current_density)
    plt.text(s=label, x=0.65 * fix_c, y=max(current_density) * 0.6)
    plt.title('scaled ref: fixed_{}={}, BAD={:.1f}, 2 params'.format(fixed_allele, fix_c, BAD))
    plt.savefig(os.path.expanduser(
        '~/fixed_alt/abcd/scaled_2params_q15-q95_{}_BAD={:.1f}_fixed_{}.png'.format(fixed_allele, BAD, fix_c)))
    plt.close(fig)
