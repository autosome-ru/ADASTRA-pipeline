import os
import sys
import numpy as np
from scipy import stats as st
from scipy import optimize
import pandas as pd
import subprocess

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import parameters_path
from scripts.HELPERS.helpers import states


def make_binom_matrix(size_of_counts, nonzero_dict, p, filename):
    if os.path.isfile(filename + '_binom.precalc.npy'):
        binom = np.load(filename + '_binom.precalc.npy')
    else:
        print('n_max = {}'.format(size_of_counts))
        binom = np.zeros((size_of_counts, size_of_counts), dtype=np.float128)
        for n in range(10, size_of_counts):
            if n not in nonzero_dict:
                continue
            print(n)
            if p != 0.5:
                f1 = st.binom(n, p).pmf
                f2 = st.binom(n, 1 - p).pmf
                binom_norm = 1 - sum(0.5 * (f1(k) + f2(k)) for k in [0, 1, 2, 3, 4, n - 4, n - 3, n - 2, n - 1, n])
                for k in range(5, n - 4):
                    binom[n, k] = 0.5 * (f1(k) + f2(k)) / binom_norm
            else:
                f = st.binom(n, p).pmf
                binom_norm = 1 - sum(f(k) for k in [0, 1, 2, 3, 4, n - 4, n - 3, n - 2, n - 1, n])
                for k in range(5, n - 4):
                    binom[n, k] = f(k) / binom_norm
        np.save(filename + '_binom.precalc.npy', binom)

    binom_sum = binom.cumsum(axis=1)
    for n in range(size_of_counts):
        binom_sum[n, n - 5:] = 1
    np.save(filename + '_binom_sum.precalc.npy', binom_sum)

    if os.path.isfile(filename + '.precalc.npy'):
        noise = np.load(filename + '.precalc.npy')
    else:
        noise = np.zeros((size_of_counts, size_of_counts), dtype=np.float128)
        for n in range(10, size_of_counts):
            if n not in nonzero_dict:
                continue
            print(n)
            for k in range(5, n - 4):
                noise[n, k] = get_noise_density(n, k)
        np.save(filename + '.precalc.npy', noise)

    noise_sum_alt = noise.cumsum(axis=1)
    for n in range(size_of_counts):
        noise_sum_alt[n, n - 5:] = 1
    np.save(filename + '_binom_sum.precalc.npy', noise_sum_alt)

    noise_sum_ref = noise.cumsum(axis=1)
    for n in range(size_of_counts):
        noise_sum_ref[n, n - 5:] = 1
    np.save(filename + '_binom_sum.precalc.npy', noise_sum_ref)

    return binom, noise


def make_counts_matrix_and_nonzero_dict(stats_pandas_dataframe):
    max_cover_in_stats = max(stats['cover'])
    counts_matrix = np.zeros((max_cover_in_stats + 1, max_cover_in_stats + 1), dtype=np.int64)
    nonzero_dict = {}

    for index, row in stats_pandas_dataframe.iterrows():
        n, k, SNP_counts = row['cover'], row['ref_counts'], row['counts']
        if k <= 4 or k >= n - 4:
            continue
        try:
            nonzero_dict[n].append(k)
        except KeyError:
            nonzero_dict[n] = [k]
        counts_matrix[n, k] = SNP_counts
    return counts_matrix, nonzero_dict


def get_sensible_n_array(n_array, samples):
    return [n for n in n_array if samples[n] >= n + 1]


def get_noise_density(n, k):
    if n == 10:
        return 1
    if k <= 4 or k >= n - 4:
        return 0
    return max(k - n / 2, 0) / (1 / 2 * (n - 5 - n // 2) * (n // 2 - 4))


def make_derivative_nonzero(counts_matrix_row, binom_matrix_row, noise_matrix_row):
    def target(alpha):
        return -1 * sum(
            sum(counts_matrix_row[k] * (
                    binom_matrix_row[k] * (-1) +
                    noise_matrix_row[k])
                / (binom_matrix_row[k] * (1 - alpha) +
                   noise_matrix_row[k] * alpha)
                for k in range(len(counts_matrix_row)) if counts_matrix_row[k] > 0))

    return target


def fit_alpha(noise_matrix_row, binom_matrix_row, counts_matrix_row):
    try:
        alpha_coefficient = optimize.brenth(
            f=make_derivative_nonzero(counts_matrix_row, binom_matrix_row, noise_matrix_row),
            a=0.01, b=0.999)
        if alpha_coefficient <= 0.001:
            return 'NaN'
    except ValueError:
        return 'NaN'
    return alpha_coefficient


def fit_weights_for_n_array(n_array, counts_matrix, nonzero_dict, BAD):
    binom_matrix, noise_matrix = make_binom_matrix(counts_matrix.shape[0], nonzero_dict, 1/(BAD + 1), filename)
    weights_of_correction = {}

    #for n in n_array:
    #    weights_of_correction[n] = fit_alpha(counts_matrix_row=counts_matrix[n, :],
    #                                         binom_matrix_row=binom_matrix[n, :],
    #                                         noise_matrix_row=noise_matrix[n, :])
    #return weights_of_correction, binom_matrix, noise_matrix
    return None, None


def make_r_lowess(weights_dict, n_array, BAD, samples, span):
    sv = pd.DataFrame({'alpha': [weights_dict[x] for x in weights_dict],
                       'snps': np.log(samples[n_array])})
    sv.index = n_array
    sv.to_csv(parameters_path + 'weights_BAD={:.1f}.tsv'.format(BAD), sep='\t')

    bash_command = "Rscript {} {:.1f} {}".format('betabin_fit.R', BAD, str(span))
    process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE)
    process.communicate()

    sv_r = pd.read_table(os.path.expanduser('~/r_weights_BAD={:.1f}.tsv'.format(BAD)))
    sv_r.index = n_array
    sv_r['lowess_r_weights'] = sv_r['lowess_r_weights'].apply(lambda x: max(x, 0))

    return dict(zip(n_array, sv_r['lowess_r_weights']))


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
        elif n > n_max:
            rv = last_w
        else:
            ns_to_fix.append(n)
            continue
        if ns_to_fix:
            for n_fix in ns_to_fix:
                extrapolated_weights[n_fix] = 0.5 * (last_w + rv)
            ns_to_fix = []
        last_w = rv
        extrapolated_weights[n] = rv
    return np.array(extrapolated_weights)


if __name__ == '__main__':
    for BAD in states:
        filename = parameters_path + 'cover_bias_statistics_BAD={:.1f}.tsv'.format(BAD)
        stats = pd.read_table(filename)
        stats['cover'] = stats['cover'].astype(int)
        stats['ref_counts'] = stats['ref_counts'].astype(int)
        counts, dict_of_nonzero_N = make_counts_matrix_and_nonzero_dict(stats)

        max_n = max(stats['cover'])
        observations = counts.sum(axis=1)
        sensible_n_array = get_sensible_n_array(range(10, max_n + 1), observations)
        weights, binom, noise = fit_weights_for_n_array(sensible_n_array, counts, dict_of_nonzero_N, BAD)

        # weights = make_r_lowess(weights, [n for n in sensible_n_array if weights[n] != 'NaN'], BAD, observations, 0.15)
        # weights = extrapolate_weights(weights, max_n)
        # np.save(parameters_path + 'weights_BAD={:.1f}.npy'.format(BAD), weights)

