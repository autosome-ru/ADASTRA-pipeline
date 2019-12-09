import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from scipy import optimize
from scipy import stats as st
import itertools


def make_binom_matrix(valid_n, p):
    n_max_from_valid_n = max(valid_n)
    rv = np.zeros((n_max_from_valid_n + 1, n_max_from_valid_n + 1), dtype=np.float128)
    noise = np.zeros((n_max_from_valid_n + 1, n_max_from_valid_n + 1), dtype=np.float128)
    for n in valid_n:
        print(n)
        # FIXME ТАК ХОТЕЛ САНЯ, ПОДУМОТЬ КАК ЛУЧЧШЕШЕШЕШЕ
        if p != 0.5:
            f1 = st.binom(n, p).pmf
            f2 = st.binom(n, 1 - p).pmf
        else:
            # f1 = st.binom(n, p).pmf
            # f2 = f1
            f = st.binom(n, p).pmf
        for k in range(n + 1):
            if p != 0.5:
                rv[n, k] = 0.5 * (f1(k) + f2(k))
                noise[n, k] = 2 * k / (n * (n + 1)) if n != 0 else 0
            else:
                rv[n, k] = f(k)
    return rv, noise


def make_counts_matrix_and_nonzero_dict(stats_pandas_dataframe):
    max_cover_in_stats = max(stats['cover'])
    counts_matrix = np.zeros((max_cover_in_stats + 1, max_cover_in_stats + 1), dtype=np.int64)
    nonzero_dict = {}

    for index, row in stats_pandas_dataframe.iterrows():
        n, k, SNP_counts = row['cover'], row['ref_counts'], row['counts']
        try:
            nonzero_dict[n].append(k)
        except KeyError:
            nonzero_dict[n] = [k]
        counts_matrix[n, k] = SNP_counts
    return counts_matrix, nonzero_dict


def make_derivative_nonzero(counts_matrix, binom_matrix, noise_matrix, window):
    #TODO: DOMASHNEE ZADANIE
    """
    :param counts_matrix:
    :param binom_matrix:
    :param noise_matrix:
    :param window:
    :return:
    """
    def target(alpha):
        return -1 * sum(
            sum(counts_matrix[n, k] * (binom_matrix[n, k] * (-1) + noise_matrix[n, k]) /
                (binom_matrix[n, k] * (1 - alpha) + noise_matrix[n, k] * alpha)
                for k in window[n])
            for n in window)

    return target


def fit_alpha(noise_matrix, binom_matrix, counts_matrix, window):
    try:
        alpha_coefficient = optimize.brenth(
            f=make_derivative_nonzero(counts_matrix, binom_matrix, noise_matrix, window), a=0.01, b=0.999)
    except ValueError:
        f = make_derivative_nonzero(counts_matrix, binom_matrix, noise_matrix, window)
        if np.sign(f(0)) == np.sign(f(0.9999)):
            if np.sign(f(0.9999)) > 0:
                return 1
            if np.sign(f(0)) < 0:
                return 0
            else:
                return 'NaN'
        else:
            return 1
    return alpha_coefficient


def fit_weights_for_n_array(n_array, counts, nonzero_dict, BAD):
    print('made ncr and noize')
    valid_n = make_array_for_valid_N(n_array, nonzero_dict)
    binom_matrix, noise = make_binom_matrix(valid_n, 1 / (BAD + 1))
    for n in valid_n:
        weights[n] = fit_alpha(counts_matrix=counts,
                               binom_matrix=binom_matrix,
                               noise_matrix=noise, window=get_window(n, nonzero_dict))
    # print(weights)
    return weights


def make_array_for_valid_N(n_array, nonzero_dict):
    valid_n = set()
    for n in n_array:
        valid_n |= get_window(n, nonzero_dict).keys()
    return valid_n


def get_window(n, nonzero_dict):
    if mode == "up_window":
        return get_window_up(n, nonzero_dict)
    else:
        return


def get_window_up(n, nonzero_dict):
    window = {}
    for key in nonzero_dict:
        if n <= key:
            window[key] = nonzero_dict[key]
    return window


if __name__ == '__main__':
    n_min = 16
    n_max = 301
    BAD = 2
    mode = "up_window"
    stats = pd.read_table('~/cover_bias_statistics_triploids.tsv')
    stats['cover'] = stats['cover'].astype(int)
    stats['ref_counts'] = stats['ref_counts'].astype(int)
    counts, dict_of_nonzero_N = make_counts_matrix_and_nonzero_dict(stats)
    print('made counts')


    # s_ns = list(range(n_min, min(n_max, 500))) + (list(range(500, 5000, 50)) if n_max > 500 else [])
    # s_ns = (list(range(n_min, 501, 10)) if n_min <= 500 else []) + list(range(max(501, n_min), n_max, 5))
    # s_ns = list(range(1000, 1501, 100))
    # s_ns = list(range(n_min, 300))
    s_ns = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200]

    plot_histograms = True
    count_stat_significance = False
    calculate_weights = False

    if calculate_weights:
        weights = fit_weights_for_n_array(s_ns, counts, dict_of_nonzero_N, BAD)

        # plot_ns = [n for n in true_ns if n <= 300]
        #
        # plt.scatter(plot_ns, [weights[k] for k in plot_ns])
        # plt.grid(True)
        # plt.xlabel('cover')
        # plt.ylabel('weight of correction')
        # # plt.title('Weight of correction ML fit on BAD=1\ncurated diploids, window +- {}'.format(window))
        # # plt.title('Weight of correction ML fit on BAD=1\ncurated diploids, cumulative up to {}'.format(n_max - 1))
        # plt.title('Weight of correction ML fit on BAD={}}\nall datasets, window +- {}'.format(BAD, window))
        # plt.show()

    # if plot_histograms:
    #     for n in s_ns:
    #         # statsplot = pd.DataFrame(stats.loc[stats['cover'] == n])
    #         # statsplot['ref_counts'] = statsplot['ref_counts'].astype(int)
    #         # statsplot['counts'] = statsplot['counts'].astype(int)
    #         # statsplot['color'] = np.abs(statsplot['ref_counts'] - n / 2)
    #         print('made data for n={}'.format(n))
    #         total_snps = sum(counts[n, k] for k in range(n + 1))
    #
    #         fig, ax = plt.subplots(figsize=(10, 8))
    #         sns.barplot(x=list(range(n + 1)), y=counts[n, 0:n + 1] / total_snps, ax=ax)
    #         plt.axvline(x=n / 2, color='black')
    #
    #         if calculate_weights:
    #             w = weights[n]
    #             print(w)
    #         else:
    #             w = 0.4
    #
    #         if mode == 'alpha':
    #             # norm = sum(st.binom(n, 0.5).pmf(x) * (1 - w) + 2 * x / (n * (n + 1)) * w for x in range(3, n - 2))
    #             # density = [0] * 3 + \
    #             #           [(st.binom(n, 0.5).pmf(x) * (1 - w) + 2 * x / (n * (n + 1)) * w) / norm for x in
    #             #            range(3, n - 2)] + \
    #             #           [0] * 3
    #             # label = 'weight of linear noize: {}\ntotal observations: {}'.format(round(w, 2), total_snps)
    #             norm = sum((0.5 * (st.binom(n, 0.33).pmf(x) + st.binom(n, 0.67).pmf(x)) * (1 - w) +
    #                         2 * x / (n * (n + 1)) * w) for x in range(3, n - 2))
    #             density = [0] * 3 + \
    #                       [(0.5 * (st.binom(n, 0.33).pmf(x) + st.binom(n, 0.67).pmf(x)) * (1 - w) + 2 * x / (
    #                               n * (n + 1)) * w) / norm for x in
    #                        range(3, n - 2)] + [0] * 3
    #             label = 'total observations: {}'.format(total_snps)
    #         elif mode == 'p':
    #             norm = sum(st.binom(n, w).pmf(x) for x in range(3, n - 2))
    #             density = [0] * 3 + [(st.binom(n, w).pmf(x)) / norm for x in range(3, n - 2)] + [0] * 3
    #             label = 'fitted binomial p: {}\ntotal observations: {}'.format(round(w, 2), total_snps)
    #
    #         plt.plot(list(range(n + 1)), density)
    #         plt.text(s=label, x=0.65 * n, y=max(density) * 0.6)
    #
    #         plt.title('ref-alt bias for BAD=2 n={}'.format(n))
    #         ax.legend().remove()
    #         plt.ylabel('count')
    #         plt.xlabel('ref_read_counts')
    #         plt.savefig(os.path.expanduser('~/ref-alt_bias_BAD=2_w=04_cn-{}.png'.format(n)))
    #         # plt.show()

    if count_stat_significance:
        p_values = []

        for ind, n in enumerate(nrange):
            w = weights[ind]
            v_list = list(itertools.chain.from_iterable([k] * int(counts[n, k]) for k in range(n + 1)))
            if mode == 'p':
                p_values.append(
                    -np.log10(st.kstest(v_list,
                                        lambda x: st.binom(n, w).cdf(x))[1]))
            elif mode == 'alpha':
                p_values.append(
                    -np.log10(st.kstest(v_list,
                                        lambda x: st.binom(n, 0.5).cdf(x) * (1 - w) + x * (x + 1) / (n * (n + 1)) * w)[
                                  1]))

        print(p_values)

        plt.scatter(nrange, p_values)
        plt.grid(True)
        plt.xlabel('cover')
        plt.ylabel('goodness of fit (-log10 p-value)')
        if mode == 'p':
            plt.title('KS-test of Binomial p ML fit on BAD=1')
        elif mode == 'alpha':
            plt.title('KS-test of linear density correction fit on BAD=1')

        # plt.show()

        # print(optimize.minimize(fun=make_target(counts, nck, n_min, n_max), x0=0, method='TNC',
        #                        jac=make_derivative(counts, nck, n_min, n_max), bounds=[(0, 0.5)]))
