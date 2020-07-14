import os
import sys
import numpy as np
import pandas as pd
#from matplotlib import pyplot as plt, ticker
from scipy import optimize
from scipy import stats as st
#import seaborn as sns

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths_for_components import parameters_path
from scripts.HELPERS.helpers import states

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


def make_counts_array_and_nonzero_set(stats_pandas_dataframe):
    max_cover_in_stats = max(stats_pandas_dataframe['{}_counts'.format(main_allele)])
    counts_array = np.zeros(max_cover_in_stats + 1, dtype=np.int64)
    nonzero_set = set()

    for index, row in stats_pandas_dataframe.iterrows():
        k, SNP_counts = row['{}_counts'.format(main_allele)], row['counts']
        nonzero_set.add(k)

        counts_array[k] = SNP_counts
    return counts_array, nonzero_set


def make_scaled_counts(stats_pandas_dataframe):
    try:
        max_cover_in_stats = max(stats_pandas_dataframe['{}_counts'.format(main_allele)])
    except ValueError as e:
        return [], set()
    counts_array = np.zeros(max_cover_in_stats + 1, dtype=np.int64)
    nonzero_set = set()

    for index, row in stats_pandas_dataframe.iterrows():
        k, SNP_counts = row['{}_counts'.format(main_allele)], row['counts']

        # if main_allele == 'ref':
        #     if k <= 4:
        #         continue
        #     k_raw = scaling[k]
        #     k_floor = np.floor(k_raw)
        #     k_ceil = np.ceil(k_raw)
        #     part = k_raw - k_floor
        #     floor_counts = np.ceil(SNP_counts * (1 - part))
        #     ceil_counts = SNP_counts - floor_counts
        #     nonzero_set.add(k_ceil)
        #     nonzero_set.add(k_floor)
        #
        #     print(k_raw, k_floor, k_ceil, SNP_counts, floor_counts, ceil_counts, part)
        #
        #     # counts_array[int(k_floor)] += SNP_counts
        #     counts_array[int(k_floor)] += floor_counts
        #     counts_array[int(k_ceil)] += ceil_counts
        # else:
        nonzero_set.add(k)
        counts_array[k] += SNP_counts
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
        cdf = lambda x: w * st.nbinom(r, get_p()).cdf(x) + (1 - w) * st.nbinom(fix_c, 1 - get_p()).cdf(x)
        q_5_l = [x for x in range(number + 1) if cdf(x) >= 1 - 1 / (10 ** 5)]
        q_5 = min(q_5_l) if q_5_l else False

        asb = counts_array[q_5:].sum() if q_5 else 0
        gof = calculate_gof(counts_array, w, r)

        label = 'negative binom fit for {}\ntotal observations: {}\nr={:.2f}, p={:.2f}, q={}, w={:.2f}\n' \
                'asb_snps={}\nasb/total={}%\ngof={}'.format(
            main_allele, total_snps, r, get_p(), q_left, w, asb, round(asb/total_snps*100, 4), round(gof, 4))
        plt.plot(list(range(n + 1)), current_density)
        plt.text(s=label, x=0.65 * n, y=max(current_density) * 0.6)
        plt.axvline(x=left_most, c='black', linestyle='--')
        plt.axvline(x=min(right_most, n), c='black', linestyle='--')
    plt.title('scaled ref: fixed_{}={}, BAD={:.1f}, 2 params'.format(fixed_allele, fix_c, BAD))
    plt.savefig(os.path.expanduser(
        '~/fixed_alt/abcd/scaled_2params_q15-q95_{}_BAD={:.1f}_fixed_{}.png'.format(fixed_allele, BAD, fix_c)))
    plt.close(fig)


def fit_negative_binom(n, counts_array):
    # f = make_derivative_nonzero(counts_matrix, binom_matrix, noise_matrix, window)
    # plot_target_function(f)
    #print("fit started")
    try:
        x = optimize.minimize(fun=make_log_likelihood(n, counts_array),
                              #x0=np.array([(fix_c if main_allele == 'alt' else scaling[fix_c], 0.5)]),
                              x0=np.array([fix_c, 0.5]),
                              bounds=[(0.00001, None), (0, 1)])
    except ValueError:
        return 'NaN', 10
    r, w = x.x
    gof = calculate_gof(counts_array, w, r)
    print(q_left, gof)
    return x, gof


def make_log_likelihood(n, counts_array):
    #print("target_made")

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


def get_p():
    return 1 / (BAD + 1)


def calculate_gof(counts_array, w, r):
    observed = counts_array.copy()
    observed[:q_left] = 0
    norm = observed.sum()
    expected = make_negative_binom_density(r, get_p(), w, len(observed) - 1) * norm

    idxs = (observed != 0) & (expected != 0)
    df = idxs.sum() - 3

    stat = np.sum(observed[idxs] * np.log(observed[idxs] / expected[idxs])) * 2

    if norm <= 1:
        return None
    else:
        score = np.sqrt(max(stat - df, 0) / (df * (norm - 1)))

    return score


if __name__ == '__main__':
    for main_allele in ("alt", "ref"):
        fixed_allele = "ref" if main_allele == "alt" else "alt"
        alleles = ('ref', 'alt')

        fix_c_array = list(range(5, 501))

        for BAD in states:

            save_array = np.zeros((max(fix_c_array) + 1, 4), dtype=np.float_)

            #scale_df = pd.read_table(os.path.expanduser('~/ref_counts_scaling_BAD={:.1f}.tsv'.format(BAD)))
            #scaling = dict(zip(scale_df['allele_reads'], scale_df['new_allele_reads']))

            # filename = os.path.expanduser('~/cover_bias_statistics_norm_diploids.tsv'.format(BAD))
            filename = parameters_path + 'fixed_alt_bias_statistics_BAD={:.1f}.tsv'.format(BAD)
            stats = pd.read_table(filename)
            for allele in alleles:
                stats['{}_counts'.format(allele)] = stats['{}_counts'.format(allele)].astype(int)

            for fix_c in fix_c_array:
                stats_filtered = stats[stats['{}_counts'.format(fixed_allele)] == fix_c]
                counts, set_of_nonzero_n = make_scaled_counts(stats_filtered)
                if len(set_of_nonzero_n) == 0 or counts.sum() < max(set_of_nonzero_n) - 5:
                    continue
                print('made counts')
                print('Fix {}={}'.format(fixed_allele, fix_c))
                number = len(counts) - 1

                left_most = 5
                q_left = 5
                right_most = len(counts) - 1

                calculate_negative_binom = True
                # weights = (fix_c, 0.5)
                # plot_histogram(number, counts, plot_fit=weights)
                if calculate_negative_binom:
                    weights, gof = fit_negative_binom(right_most, counts)
                    # plot_histogram(number, counts, plot_fit=weights.x)
                    save_array[fix_c, :2] = weights.x
                    save_array[fix_c, 2] = weights.success
                    save_array[fix_c, 3] = gof

            np.save(os.path.expanduser(parameters_path + 'NBweights_{}_BAD={:.1f}'.format(fixed_allele, BAD)), save_array)
