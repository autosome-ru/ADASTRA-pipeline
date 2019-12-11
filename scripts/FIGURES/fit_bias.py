import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy import optimize
from scipy import stats as st
from sklearn import metrics
import seaborn as sns


def make_binom_matrix(size_of_counts, nonzero_dict, p):
    if os.path.isfile(filename + '_binom.precalc.npy') and os.path.isfile(filename + '_linear.precalc.npy'):
        rv = np.load(filename + '_binom.precalc.npy')
        noise = np.load(filename + '_linear.precalc.npy')
        return rv, noise

    print('n_max = {}'.format(size_of_counts))
    rv = np.zeros((size_of_counts, size_of_counts), dtype=np.float128)
    noise = np.zeros((size_of_counts, size_of_counts), dtype=np.float128)
    for n in range(6, size_of_counts):
        if n not in nonzero_dict:
            continue
        print(n)
        if p != 0.5:
            f1 = st.binom(n, p).pmf
            f2 = st.binom(n, 1 - p).pmf
            binom_norm = 1 - sum(0.5 * (f1(k) + f2(k)) for k in [0, 1, 2, n - 2, n - 1, n])
            for k in range(3, n - 2):
                rv[n, k] = 0.5 * (f1(k) + f2(k)) / binom_norm
                noise[n, k] = get_noise_density(n, k)
        else:
            f = st.binom(n, p).pmf
            binom_norm = 1 - sum(f(k) for k in [0, 1, 2, n - 2, n - 1, n])
            for k in range(3, n - 2):
                rv[n, k] = f(k) / binom_norm
                noise[n, k] = get_noise_density(n, k)
    np.save(filename + '_binom.precalc.npy', rv)
    np.save(filename + '_linear.precalc.npy', noise)
    return rv, noise


def get_noise_density(n, k):
    return 2 * k / (n * (n - 5)) if n != 0 else 0


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
    """
    :param counts_matrix: number of counts 2D np.array
    :param binom_matrix:
    :param noise_matrix:
    :param window: dict[n] = array of valid k
    :return: target function for optimization
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
            if np.sign(f(0.9999)) < 0:
                return 1
            if np.sign(f(0)) > 0:
                return 0
            else:
                return 'NaN'
        else:
            return 1
    return alpha_coefficient


def plot_target_function(f):
    x = [v / 100 for v in range(1, 100)]
    y = [f(v) for v in x]
    plt.scatter(x, y)
    plt.grid(True)
    plt.show()


def fit_weights_for_n_array(n_array, counts_matrix, nonzero_dict, samples):
    print('made ncr and noize')
    binom_matrix, noise = make_binom_matrix(counts_matrix.shape[0], nonzero_dict, get_p())
    print(binom_matrix[15, :16])
    weights_of_correction = {}
    for n in n_array:
        print('fitting for n={}'.format(n))
        weights_of_correction[n] = fit_alpha(counts_matrix=counts_matrix,
                                             binom_matrix=binom_matrix,
                                             noise_matrix=noise, window=get_window(n, nonzero_dict, samples))
        print(weights_of_correction[n])
    return weights_of_correction


def get_window(n, nonzero_dict, samples, window_mode=None):
    if window_mode is None:
        window_mode = mode
    if window_mode == "up_window":
        return get_window_up(n, nonzero_dict)
    elif window_mode == "up_window_n_sq":
        return get_window_up_n_sq(n, nonzero_dict, samples)
    else:
        return


def get_window_up(n, nonzero_dict):
    window = {}
    for key in nonzero_dict:
        if n <= key:
            window[key] = nonzero_dict[key]
    return window


def get_window_up_n_sq(n, nonzero_dict, samples):
    window = {}
    current_cumulative_counts = 0
    required_counts = (n + 1) ** 2
    for key in sorted(list(nonzero_dict.keys())):
        if n <= key:
            window[key] = nonzero_dict[key]
            current_cumulative_counts += samples[key]
        if current_cumulative_counts >= required_counts:
            break
    return window


def plot_fit(weights_of_correction, save=True):
    fig, ax = plt.subplots(figsize=(10, 8))
    plt.scatter(list(weights_of_correction.keys()), [weights_of_correction[k] for k in weights_of_correction])
    plt.grid(True)
    plt.xlabel('cover')
    plt.ylabel('weight of correction')
    plt.title('Weight of correction ML fit on BAD={}\nall_datasets, {}'.format(BAD, mode))
    if save:
        plt.savefig(os.path.expanduser('~/plots/weights_BAD={}_mode={}.png'.format(BAD, mode)))
    else:
        plt.show()


def plot_quality(scores, binom_scores, metric_mode, save=True):
    print(scores, binom_scores)
    fig, ax = plt.subplots(figsize=(10, 8))
    plt.scatter(list(scores.keys()), [scores[k] for k in scores], label='fit')
    plt.scatter(list(binom_scores.keys()), [binom_scores[k] for k in binom_scores], label='binom')
    max_y = max([binom_scores[k] for k in binom_scores])
    ax.set_ylim([-1 * max_y / 50, max_y * 1.02])
    plt.grid(True)
    plt.xlabel('cover')
    plt.ylabel('Scores')
    plt.title('Scores of correction ML fit on BAD={}\nall_datasets, {}'.format(BAD, mode))
    plt.legend()
    if save:
        plt.savefig(os.path.expanduser('~/plots/qual_of_fit_{}_BAD={}_mode={}.png'.format(metric_mode, BAD, mode)))
    else:
        plt.show()


def calculate_score(weights_of_correction, counts_matrix, metric_mode):
    scores = dict()
    binom_scores = dict()
    for n in weights_of_correction:
        norm, observed = get_observed(n, counts_matrix, normalize=False)
        if norm == 0:
            continue
        expected = get_probability_density(n, weights_of_correction[n]) * norm
        expected_binom = get_probability_density(n, 0) * norm

        print(n, norm, norm >= (n + 1) ** 2)

        # plot_distributions(n, expected, expected_binom, observed)

        if metric_mode == 'rmse':
            scores[n] = np.sqrt(metrics.mean_squared_error(expected, observed))
            binom_scores[n] = np.sqrt(metrics.mean_squared_error(expected_binom, observed))
        elif metric_mode == 'chi_sq':
            scores[n] = st.chisquare(observed[3: n - 2], expected[3: n - 2])[1]
            binom_scores[n] = st.chisquare(observed[3: n - 2], expected_binom[3: n - 2])[1]
        elif metric_mode == 'g':
            scores[n] = np.sum(observed[3: n - 2] * np.log(observed[3: n - 2] / expected[3: n - 2])) * 2
            binom_scores[n] = np.sum(observed[3: n - 2] * np.log(observed[3: n - 2] / expected_binom[3: n - 2])) * 2

    return scores, binom_scores


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


def get_probability_density(n, alpha):
    p = get_p()
    f1 = st.binom(n, p).pmf
    f2 = st.binom(n, 1 - p).pmf
    norm = sum((0.5 * (f1(x) + f2(x)) * (1 - alpha) +
                2 * x / (n * (n + 1)) * alpha) for x in range(3, n - 2))
    density = [0] * 3 + \
              [(0.5 * (f1(x) + f2(x)) * (1 - alpha) +
                2 * x / (n * (n + 1)) * alpha) / norm for x in range(3, n - 2)] + [0] * 3
    return np.array(density)


def plot_window_sizes_in_snps(n_array, nonzero_dict, samples, window_mode, save=True):
    y = [np.log10(sum(samples[n] for n in get_window(n, nonzero_dict, samples, window_mode))) for n in n_array]
    fig, ax = plt.subplots(figsize=(10, 8))
    plt.scatter(n_array, y, label='in window')
    plt.scatter(n_array, np.log10(np.array(n_array) + 1), s=1, label='linear')
    plt.scatter(n_array, 2 * np.log10(np.array(n_array) + 1), s=1, label='sq')
    plt.grid(True)
    plt.xlabel('cover')
    plt.ylabel('log10 number of snps in window')
    plt.title('Number of observations in window on BAD={}\nall_datasets, {}'.format(BAD, window_mode))
    plt.legend()
    if save:
        plt.savefig(os.path.expanduser('~/plots/obs_in_window_BAD={}_mode={}.png'.format(BAD, window_mode)))
    else:
        plt.show()


def plot_histogram(n, weight, save=True):
    print('made data for n={}'.format(n))
    current_density = get_probability_density(n, weight)
    total_snps = counts[n, 0:n + 1].sum()

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.barplot(x=list(range(n + 1)), y=counts[n, 0:n + 1] / total_snps, ax=ax)
    plt.axvline(x=n / 2, color='black')

    label = 'weight of linear noize: {}\ntotal observations: {}'.format(round(weight, 2),
                                                                        total_snps)
    plt.plot(list(range(n + 1)), current_density)
    plt.text(s=label, x=0.65 * n, y=max(current_density) * 0.6)

    plt.title('ref-alt bias for BAD={} n={}'.format(BAD, n))
    ax.legend().remove()
    plt.ylabel('count')
    plt.xlabel('ref_read_counts')
    if save:
        plt.savefig(os.path.expanduser('~/plots/ref-alt_bias_n={}_BAD={}_w={}.png'.format(n, BAD, weight)))
    else:
        plt.show()


def get_max_sensible_n(n_array, samples, nonzero_dict, sensible_mode='linear'):
    for i, n in enumerate(n_array):
        current_cumulative_counts = 0
        if sensible_mode == 'linear':
            required_counts = (n + 1)
        elif sensible_mode == 'sq':
            required_counts = (n + 1) ** 2
        else:
            raise ValueError(sensible_mode)

        for key in sorted(list(nonzero_dict.keys())):
            if n <= key:
                current_cumulative_counts += samples[key]
            if current_cumulative_counts >= required_counts:
                break
        else:
            return n_array[max(0, i - 1)]
    return n_array[-1]


if __name__ == '__main__':

    BAD = 1
    mode = "up_window_n_sq"
    metric_modes = ['rmse', 'g']
    filename = os.path.expanduser('~/cover_bias_statistics_norm_diploids.tsv')
    stats = pd.read_table(filename)
    stats['cover'] = stats['cover'].astype(int)
    stats['ref_counts'] = stats['ref_counts'].astype(int)
    counts, dict_of_nonzero_N = make_counts_matrix_and_nonzero_dict(stats)
    total_snps_with_cover_n = counts.sum(axis=1)
    print('made counts')

    # s_ns = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200]
    # s_ns = range(10, max(stats['cover']) + 1, 5)
    s_ns = range(6, 231, 1)

    max_sensible_n = get_max_sensible_n(s_ns, total_snps_with_cover_n, dict_of_nonzero_N)

    plot_window_counts = True

    calculate_weights = True
    plot_fit_weights = True

    calculate_fit_quality = True
    plot_fit_quality = True
    
    plot_histograms = True

    if plot_window_counts:
        for window_mode in ("up_window", "up_window_n_sq"):
            plot_window_sizes_in_snps(s_ns, dict_of_nonzero_N, total_snps_with_cover_n, window_mode)

    if calculate_weights:
        sensible_n_array = [n for n in s_ns if n <= max_sensible_n]
        weights = fit_weights_for_n_array(sensible_n_array, counts, dict_of_nonzero_N, total_snps_with_cover_n)
        if plot_fit_weights:
            plot_fit(weights)

        if calculate_fit_quality:
            for metric in metric_modes:
                calculated_fit_metrics, calculated_binom_metrics = calculate_score(weights, counts, metric)

                if plot_fit_quality:
                    plot_quality(calculated_fit_metrics, calculated_binom_metrics, metric)

        if plot_histograms:
            for n in [min(sensible_n_array), get_max_sensible_n(s_ns, total_snps_with_cover_n, dict_of_nonzero_N,
                                                                'sq'), max_sensible_n]:
                plot_histogram(n, weights[n], save=True)
