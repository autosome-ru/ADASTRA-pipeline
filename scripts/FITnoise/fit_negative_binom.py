import os
import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import optimize
from scipy import stats as st
from scipy.special import beta, psi, poch, comb
from sklearn import metrics
import seaborn as sns
import subprocess

pd.set_option('display.max_columns', 7)


def make_negative_binom_density(r, p, size_of_counts):
    negative_binom_density_array = np.zeros(size_of_counts + 1, dtype=np.float128)
    f = st.nbinom(r, p).pmf
    negative_binom_norm = 1.0 - np.array([f(k) for k in [0, 1, 2, 3, 4]]).sum()
    for k in range(5, size_of_counts + 1):
        negative_binom_density_array[k] = f(k) / np.float(negative_binom_norm)
    print("negative_binom_density counted")
    return negative_binom_density_array


def make_counts_array_and_nonzero_set(stats_pandas_dataframe):
    max_cover_in_stats = max(stats_pandas_dataframe['{}_read_counts'.format(alt)])
    counts_array = np.zeros(max_cover_in_stats + 1, dtype=np.int64)
    nonzero_set = set()

    for index, row in stats_pandas_dataframe.iterrows():
        k, SNP_counts = row['{}_read_counts'.format(alt)], row['counts']
        nonzero_set.add(k)

        counts_array[k] = SNP_counts
    return counts_array, nonzero_set


def plot_histogram(n, counts_array, plot_fit=None, save=True):
    print('made data for n={}'.format(n))

    total_snps = counts_array[0:n + 1].sum()

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.barplot(x=list(range(n + 1)),
                y=counts_array[0:n + 1] / total_snps, ax=ax)
    if plot_fit is not None:
        r = plot_fit[0]
        p = plot_fit[1]
        current_density = make_negative_binom_density(r, p, n)
        label = 'negative binom fit\ntotal observations: {}\nr={:.5f}, p={:.5f}'.format(total_snps, r, p)
        plt.plot(list(range(n + 1)), current_density)
        plt.text(s=label, x=0.65 * n, y=max(current_density) * 0.6)
    plt.show()


def fit_negative_binom(n, counts_array):
    # f = make_derivative_nonzero(counts_matrix, binom_matrix, noise_matrix, window)
    # plot_target_function(f)
    print("fit started")
    try:
        x = optimize.minimize(fun=make_log_likelihood(counts_array, n), x0=np.array([(0.001, 0.15)]),

                              bounds=[(0.00001, None), (0.1, 0.3)])
    except ValueError:
        return 'NaN'
    return x


def make_log_likelihood(n, counts_array):
    print("target_made")

    def target(x):
        r = x[0]
        p = x[1]
        print(r, p)
        print("Counting likelihood")
        neg_bin_dens = make_negative_binom_density(r, p, len(counts_array))
        return -1 * sum(counts_array[k] * np.log(neg_bin_dens[k])
                        for k in range(5, n) if counts_array[k] != 0)

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
    alt = "ref"
    for BAD in [1]:

        # filename = os.path.expanduser('~/cover_bias_statistics_norm_diploids.tsv'.format(BAD))
        filename = os.path.expanduser('~/{}_bias_statistics_BAD={:.1f}.tsv'.format(alt, BAD))
        stats = pd.read_table(filename)
        stats['{}_read_counts'.format(alt)] = stats['{}_read_counts'.format(alt)].astype(int)
        counts, dict_of_nonzero_N = make_counts_array_and_nonzero_set(stats)
        print('made counts')
        number = 40

        calculate_negative_binom = True
        # weights = (0.00001, 0.15)
        # plot_histogram(number, counts, plot_fit=weights)
        if calculate_negative_binom:
            weights = fit_negative_binom(counts, number)
            print(weights)
            plot_histogram(number, counts, plot_fit=weights.x)
