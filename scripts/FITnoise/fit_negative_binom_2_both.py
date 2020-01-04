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

pd.set_option('display.max_columns', 7)


def make_negative_binom_density(r, p, w, size_of_counts, for_plot=False):
    negative_binom_density_array = np.zeros(size_of_counts + 1, dtype=np.float128)
    dist1 = st.nbinom(r, p)
    f1 = dist1.pmf
    cdf1 = dist1.cdf
    dist2 = st.nbinom(r, 1 - p)
    f2 = dist2.pmf
    cdf2 = dist2.cdf
    negative_binom_norm = (cdf1(size_of_counts) - cdf1(left_most - 1)) * w + \
                          (cdf2(size_of_counts) - cdf2(left_most - 1)) * (1 - w)
    plot_norm = (cdf1(size_of_counts) - cdf1(4)) * w + \
                (cdf2(size_of_counts) - cdf2(4)) * (1 - w)
    for k in range(5, size_of_counts + 1):
        if for_plot:
            negative_binom_density_array[k] = (w * f1(k) + (1 - w) * f2(k)) / plot_norm
        else:
            negative_binom_density_array[k] = (w * f1(k) + (1 - w) * f2(k)) / negative_binom_norm
    return negative_binom_density_array


def make_scaled_counts(stats_ref, stats_alt):
    scale_df = pd.read_table(os.path.expanduser('~/ref_counts_scaling_BAD={:.1f}.tsv'.format(BAD)))
    scaling = dict(zip(scale_df['allele_reads'], scale_df['new_allele_reads']))
    stats_frames = {'ref': stats_ref, 'alt': stats_alt}
    nonzero_set = set()

    max_cover_in_stats = max(
        max(stats_frames[fixed_allele]['{}_counts'.format('ref' if fixed_allele == 'alt' else 'alt')])
        for fixed_allele in ('ref', 'alt'))
    counts_array = np.zeros(max_cover_in_stats + 2, dtype=np.int64)

    for fixed_allele in ('ref', 'alt'):
        main_allele = 'ref' if fixed_allele == 'alt' else 'alt'
        for index, row in stats_frames[fixed_allele].iterrows():
            k, SNP_counts = row['{}_counts'.format(main_allele)], row['counts']
            if k <= 4:
                continue
            k_raw = scaling[k]
            k_floor = np.floor(k_raw)
            k_ceil = np.ceil(k_raw)
            part = k_raw - k_floor
            floor_counts = np.ceil(SNP_counts * (1 - part))
            ceil_counts = SNP_counts - floor_counts
            nonzero_set.add(k_ceil)
            nonzero_set.add(k_floor)

            # print(k_raw, k_floor, k_ceil, SNP_counts, floor_counts, ceil_counts, part)

            # counts_array[int(k_floor)] += SNP_counts
            counts_array[int(k_floor)] += floor_counts
            counts_array[int(k_ceil)] += ceil_counts
    return counts_array, nonzero_set


def plot_histogram(n, counts_array, plot_fit=None, save=True):
    print('made data for n={}'.format(n))

    total_snps = counts_array[0:n + 1].sum()

    fig, ax = plt.subplots(figsize=(10, 8))
    x = list(range(n + 1))
    sns.barplot(x=x,
                y=counts_array[0:n + 1] / total_snps, ax=ax)
    ax.xaxis.set_major_locator(ticker.FixedLocator(np.arange(0, len(x), 5)))
    ax.xaxis.set_major_formatter(ticker.FixedFormatter(x[::5]))
    ax.tick_params(axis="x", rotation=90)
    if plot_fit is not None:
        r = plot_fit[0]
        w = plot_fit[1]
        current_density = make_negative_binom_density(r, get_p(), w, n, for_plot=True)
        label = 'negative binom fit\ntotal observations: {}\nr={:.2f}, p={:.2f},\nq{}={}, q{}={}\nw={:.2f}'.format(
            total_snps, r, get_p(), t_left, q_left, t_right, q_right, w)
        plt.plot(list(range(n + 1)), current_density)
        plt.text(s=label, x=0.65 * n, y=max(current_density) * 0.6)
        plt.axvline(x=left_most, c='black', linestyle='--')
        plt.axvline(x=min(right_most, n), c='black', linestyle='--')
    plt.title('scaled ref: fixed_allele={}, BAD={:.1f}, 2 params'.format(fix_c, BAD))
    plt.savefig(os.path.expanduser(
        '~/fixed_alt/scaled_2params_both_q{}-q{}_BAD={:.1f}_fixed_{}.png'.format(t_left, t_right, BAD, fix_c)))
    plt.close(fig)


def fit_negative_binom(n, counts_array):
    # f = make_derivative_nonzero(counts_matrix, binom_matrix, noise_matrix, window)
    # plot_target_function(f)
    print("fit started")
    try:
        x = optimize.minimize(fun=make_log_likelihood(n, counts_array),
                              x0=np.array([(fix_c, 0.1)]),

                              bounds=[(0.00001, None), (0, 0.5)])
    except ValueError:
        return 'NaN'
    return x


def make_log_likelihood(n, counts_array):
    print("target_made")

    def target(x):
        r = x[0]
        w = x[1]
        # print(r, p)
        # print("Counting likelihood")
        neg_bin_dens = make_negative_binom_density(r, get_p(), w, right_most)
        return -1 * sum(counts_array[k] * np.log(neg_bin_dens[k])
                        for k in range(left_most, n) if counts_array[k] != 0)

    return target


def plot_target_function(f):
    x = [v / 100 for v in range(1, 100)]
    y = [f(v) for v in x]
    plt.scatter(x, y)
    plt.grid(True)
    plt.show()


def plot_distributions(n, expected, expected_binom, observed):
    plt.scatter(list(range(n + 1)), expected, label='fit')
    plt.scatter(list(range(n + 1)), expected_binom, label='binom')
    plt.scatter(list(range(n + 1)), observed, label='observed')
    plt.grid(True)
    plt.legend()
    plt.show()


def get_p():
    return 1 / (BAD + 1)


def get_observed(n, counts_matrix, normalize=True):
    norm = counts_matrix[n, :n + 1].sum()
    if norm == 0:
        return 0, counts_matrix[n, :n + 1]
    observed = counts_matrix[n, :n + 1]
    if normalize:
        observed = observed / norm
    return norm, np.array(observed)


def extrapolate_weights(weights_of_correction, n_max):
    extrapolated_weights = np.zeros(n_max + 1)
    n_min = 10
    last_w = 0
    ns_to_fix = []
    for n in range(n_max):
        if n < n_min:
            rv = 0
        elif n in weights_of_correction:
            rv = weights_of_correction[n]
        elif n > max(weights_of_correction):
            extrapolated_weights[n:] = last_w
            break
        else:
            ns_to_fix.append(n)
            continue
        if ns_to_fix:
            for n_fix in ns_to_fix:
                extrapolated_weights[n_fix] = 0.5 * (last_w + rv)
            ns_to_fix = []
        last_w = rv
        extrapolated_weights[n] = rv
    # print(extrapolated_weights)
    # plt.plot(extrapolated_weights)
    # plt.show()
    return np.array(extrapolated_weights)


if __name__ == '__main__':
    fix_c_array = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100]
    save_array = np.zeros((max(fix_c_array), 2), dtype=np.float_)
    for BAD in [6]:  # [1, 4/3, 1.5, 2, 2.5, 3, 4, 5, 6]:
        for fix_c in fix_c_array:

            # filename = os.path.expanduser('~/cover_bias_statistics_norm_diploids.tsv'.format(BAD))
            filename = os.path.expanduser('~/fixed_alt_bias_statistics_BAD={:.1f}.tsv'.format(BAD))
            stats = pd.read_table(filename)
            for allele in ('ref', 'alt'):
                stats['{}_counts'.format(allele)] = stats['{}_counts'.format(allele)].astype(int)
            for allele in ('ref', 'alt'):
                stats_fixed_ref = stats[stats['ref_counts'] == fix_c]
                stats_fixed_alt = stats[stats['alt_counts'] == fix_c]
            # counts, dict_of_nonzero_N = make_counts_array_and_nonzero_set(stats)
            counts, dict_of_nonzero_N = make_scaled_counts(stats_fixed_ref, stats_fixed_alt)
            print('made counts')
            number = 40

            number = min(len(counts) - 1, 250)

            if BAD == 6:
                t_left = 10
            else:
                t_left = 15
            t_right = 100

            # cdf = lambda x: 0.5 * (st.nbinom(fix_c, get_p()).cdf(x) + st.nbinom(fix_c, 1 - get_p()).cdf(x))
            # q_left = min(x for x in range(number + 1) if cdf(x) >= t_left / 100)
            # q_right = max(x for x in range(number + 1) if cdf(x) <= t_right / 100)

            total_snps = counts.sum()
            cdf = counts.cumsum() / total_snps
            assert cdf[-1] == 1
            q_left = min(x for x in range(len(counts)) if cdf[x] >= t_left / 100)
            q_right = max(x for x in range(len(counts)) if cdf[x] <= t_right / 100)

            print('q={} {}'.format(q_left, q_right))
            left_most = max(5, q_left)
            # right_most = min(q_right, len(counts))
            right_most = len(counts)

            q_right = min(q_right, number)

            calculate_negative_binom = True
            # weights = (fix_c, 0.5)
            # plot_histogram(number, counts, plot_fit=weights)

            r = []

            if calculate_negative_binom:
                weights = fit_negative_binom(right_most, counts)
                print(weights)
                # write to saver
                save_array[fix_c, :] = weights.x
                #plot_histogram(number, counts, plot_fit=weights.x)
    np.save(filename_for_np_array + "{}".format(BAD), save_array)
