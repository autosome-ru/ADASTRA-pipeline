import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from scipy import optimize
from scipy import stats as st
import operator as op
from functools import reduce
import itertools


def ncr(n, r):
    r = min(r, n - r)
    numer = reduce(op.mul, range(n, n - r, -1), 1)
    denom = reduce(op.mul, range(1, r + 1), 1)
    return numer / denom


def make_ncr_array(n_max, ns=None):
    rv = np.zeros((n_max + 1 + window, n_max + 1 + window), dtype=np.float128)
    for n in range(n_max + 1 + window):
        if ns:
            if n not in ns: continue
        print(n)
        # n_pow = 2 ** (-n)
        f = st.binom(n, 0.5).pmf
        for k in range(n + 1):
            rv[n, k] = f(k)
    return rv


def make_ncr_array_full(n_max, p):
    n_real_max = min(max(stats['cover']), n_max + window)
    rv = np.zeros((n_real_max + 1, n_real_max + 1), dtype=np.float128)
    noise = np.zeros((n_real_max + 1, n_real_max + 1), dtype=np.float128)
    for n in range(n_min, n_max + 1 + window):
        print(n)
        # n_pow = 2 ** (-n)
        # f = st.binom(n, 0.5).pmf
        if p != 0.5:
            f1 = st.binom(n, p).pmf
            f2 = st.binom(n, 1 - p).pmf
        else:
            f1 = st.binom(n, p).pmf
            f2 = f1
        for k in range(n + 1):
            # rv[n, k] = f(k)
            rv[n, k] = 0.5 * (f1(k) + f2(k))
            noise[n, k] = 2 * k / (n * (n + 1)) if n != 0 else 0
    return rv, noise


def make_counts_array(stats, n_max, ns=None):
    c = np.zeros((n_max + 1 + window, n_max + 1 + window), dtype=np.int64)
    for n in range(n_max + 1 + window):
        if ns:
            if n not in ns: continue
        print(n)
        for k in range(n + 1):
            slice = stats[(stats['cover'] == n) & (stats['ref_counts'] == k)]
            if not slice.empty:
                assert len(slice.index) == 1
                c[n, k] = np.int64(slice['counts'])
    return c


def make_counts_array_full(n_max, stats):
    n_real_max = min(max(stats['cover']), n_max + window)
    c = np.zeros((n_real_max + 1, n_real_max + 1), dtype=np.int64)
    for index, row in stats.iterrows():
        n, k, cts = row['cover'], row['ref_counts'], row['counts']
        if n > n_real_max: continue
        c[n, k] = cts
    return c


def make_target(counts, nck, n_min, n_max):
    def target(alpha):
        return -1 * sum(
            sum(counts[n, k] * np.log(nck[n, k] * (1 - alpha) + 2 * k / (n * (n + 1)) * alpha)
                for k in range(0, n + 1))
            for n in range(n_min, n_max + 1))

    return target


def make_derivative(counts, nck, n_min, n_max):
    def target(alpha):
        return -1 * sum(
            sum(counts[n, k] * (nck[n, k] * (-1) + noise[n, k]) /
                (nck[n, k] * (1 - alpha) + noise[n, k] * alpha)
                for k in range(0, n + 1))
            for n in range(n_min, n_max + 1))

    return target


def make_derivative_nonzero(counts, nck, noise, nonzero_dict, n_min, n_max):
    def target(alpha):
        return -1 * sum(
            sum(counts[n, k] * (nck[n, k] * (-1) + noise[n, k]) /
                (nck[n, k] * (1 - alpha) + noise[n, k] * alpha)
                for k in nonzero_dict[n])
            for n in [key for key in nonzero_dict if n_min <= key <= n_max])

    return target


def make_derivative_p(counts, n_min, n_max):
    def target(p):
        return -1 * sum(
            sum(counts[n, k] * (k - n * p) / (p * (1 - p))
                for k in range(0, n + 1))
            for n in range(n_min, n_max + 1))

    return target


def make_ns(s_ns, window):
    ns = []
    for n in s_ns:
        for k in range(n - window, n + window + 1):
            if k >= n_min and k not in ns:
                ns.append(k)
    return ns


def make_nonzero_dict(c, n_max, stats):
    n_real_max = min(max(stats['cover']), n_max + window)
    nonzero = dict()
    for n in range(n_real_max + 1):
        rv = []
        for k in range(n + 1):
            if c[n, k]:
                rv.append(k)
        if rv:
            nonzero[n] = rv
    return nonzero


def fit_alpha(n_min, n_max, noise_matrix, binom_matrix, counts_matrix):
    try:
        alpha_coefficient = optimize.brenth(
            f=make_derivative_nonzero(counts_matrix, binom_matrix, noise_matrix, nonzero_dict, n_min, n_max),
            a=0.01, b=0.999)
    except ValueError:
        f = make_derivative(counts, nck, n_min_t, n_max_t)
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


def fit_weights_for_window(n_min_input, n_max_input):
    print('made ncr and noize')
    nck, noise = make_ncr_array_full(n_max, ns)
    for n in true_ns:

        ## window definition from n
        # n_max_t = n + window
        # n_max_t = n_max + window
        # n_min_t = max(n - window, n_min)
        n_min_t = n - window
        n_max_t = n + window + 1

        if ns:
            if n not in ns:
                continue
        print(n, n_min_t, n_max_t)

        weights[n] = fit_alpha(n_min=n_min_t, n_max=n_max_t, counts_matrix=counts, binom_matrix=nck,
                               noise_matrix=noise)
    print(weights)
    return weights


def make_array_of_for_valid_N():
    pass


if __name__ == '__main__':
    n_min = 16
    n_max = 301
    window = 1
    BAD = 2

    # s_ns = list(range(n_min, min(n_max, 500))) + (list(range(500, 5000, 50)) if n_max > 500 else [])
    # s_ns = (list(range(n_min, 501, 10)) if n_min <= 500 else []) + list(range(max(501, n_min), n_max, 5))
    # s_ns = list(range(1000, 1501, 100))
    # s_ns = list(range(n_min, 300))
    s_ns = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200]
    # n-s that we wish to examine

    plot_histograms = True
    count_stat_significance = False
    calculate_weights = False

    ns = make_ns(s_ns, window)
    # n-s from s-ns +- window

    print(s_ns)
    print(ns)
    stats = pd.read_table('~/cover_bias_statistics_triploids.tsv')
    stats['cover'] = stats['cover'].astype(int)
    stats['ref_counts'] = stats['ref_counts'].astype(int)
    counts = make_counts_array_full(n_max, stats)
    print('made counts')

    nonzero_dict = make_nonzero_dict(counts, n_max, stats)

    nrange = range(n_min + 1, n_max)

    true_ns = [n for n in s_ns if [key for key in nonzero_dict if max(n - window, n_min) <= key <= n_max + window]]
    # n-s from ns that happen at least once

    weights = dict()
    if calculate_weights:
        print('made ncr and noize')
        nck, noise = make_ncr_array_full(n_max, ns)
        for n in true_ns:

            ## window definition from n
            # n_max_t = n + window
            # n_max_t = n_max + window
            # n_min_t = max(n - window, n_min)
            n_min_t = n - window
            n_max_t = n + window + 1

            if ns:
                if n not in ns:
                    continue
            print(n, n_min_t, n_max_t)

            weights[n] = fit_alpha(n_min=n_min_t, n_max=n_max_t, counts_matrix=counts, binom_matrix=nck,
                                   noise_matrix=noise)
        print(weights)

        plot_ns = [n for n in true_ns if n <= 300]

        plt.scatter(plot_ns, [weights[k] for k in plot_ns])
        plt.grid(True)
        plt.xlabel('cover')
        plt.ylabel('weight of correction')
        # plt.title('Weight of correction ML fit on BAD=1\ncurated diploids, window +- {}'.format(window))
        # plt.title('Weight of correction ML fit on BAD=1\ncurated diploids, cumulative up to {}'.format(n_max - 1))
        plt.title('Weight of correction ML fit on BAD={}}\nall datasets, window +- {}'.format(BAD, window))
        plt.show()

    if plot_histograms:
        for n in s_ns:
            # statsplot = pd.DataFrame(stats.loc[stats['cover'] == n])
            # statsplot['ref_counts'] = statsplot['ref_counts'].astype(int)
            # statsplot['counts'] = statsplot['counts'].astype(int)
            # statsplot['color'] = np.abs(statsplot['ref_counts'] - n / 2)
            print('made data for n={}'.format(n))
            total_snps = sum(counts[n, k] for k in range(n + 1))

            fig, ax = plt.subplots(figsize=(10, 8))
            sns.barplot(x=list(range(n + 1)), y=counts[n, 0:n + 1] / total_snps, ax=ax)
            plt.axvline(x=n / 2, color='black')

            if calculate_weights:
                w = weights[n]
                print(w)
            else:
                w = 0.4

            if mode == 'alpha':
                # norm = sum(st.binom(n, 0.5).pmf(x) * (1 - w) + 2 * x / (n * (n + 1)) * w for x in range(3, n - 2))
                # density = [0] * 3 + \
                #           [(st.binom(n, 0.5).pmf(x) * (1 - w) + 2 * x / (n * (n + 1)) * w) / norm for x in
                #            range(3, n - 2)] + \
                #           [0] * 3
                # label = 'weight of linear noize: {}\ntotal observations: {}'.format(round(w, 2), total_snps)
                norm = sum((0.5 * (st.binom(n, 0.33).pmf(x) + st.binom(n, 0.67).pmf(x)) * (1 - w) +
                            2 * x / (n * (n + 1)) * w) for x in range(3, n - 2))
                density = [0] * 3 + \
                          [(0.5 * (st.binom(n, 0.33).pmf(x) + st.binom(n, 0.67).pmf(x)) * (1 - w) + 2 * x / (
                                  n * (n + 1)) * w) / norm for x in
                           range(3, n - 2)] + [0] * 3
                label = 'total observations: {}'.format(total_snps)
            elif mode == 'p':
                norm = sum(st.binom(n, w).pmf(x) for x in range(3, n - 2))
                density = [0] * 3 + [(st.binom(n, w).pmf(x)) / norm for x in range(3, n - 2)] + [0] * 3
                label = 'fitted binomial p: {}\ntotal observations: {}'.format(round(w, 2), total_snps)

            plt.plot(list(range(n + 1)), density)
            plt.text(s=label, x=0.65 * n, y=max(density) * 0.6)

            plt.title('ref-alt bias for BAD=2 n={}'.format(n))
            ax.legend().remove()
            plt.ylabel('count')
            plt.xlabel('ref_read_counts')
            plt.savefig(os.path.expanduser('~/ref-alt_bias_BAD=2_w=04_cn-{}.png'.format(n)))
            # plt.show()

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
