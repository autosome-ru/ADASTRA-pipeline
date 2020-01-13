import os
import sys
import numpy as np
import pandas as pd
from scipy import optimize
from scipy import stats as st
from sklearn.linear_model import LinearRegression

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


def make_counts(stats_pandas_dataframe):
    try:
        max_cover_in_stats = max(stats_pandas_dataframe['{}_counts'.format(main_allele)])
    except ValueError as e:
        return [], set()
    counts_array = np.zeros(max_cover_in_stats + 1, dtype=np.int64)
    nonzero_set = set()

    for index, row in stats_pandas_dataframe.iterrows():
        k, SNP_counts = row['{}_counts'.format(main_allele)], row['counts']

        nonzero_set.add(k)
        counts_array[k] += SNP_counts
    return counts_array, nonzero_set


def fit_negative_binom(n, counts_array):
    # f = make_derivative_nonzero(counts_matrix, binom_matrix, noise_matrix, window)
    # plot_target_function(f)
    #print("fit started")
    try:
        x = optimize.minimize(fun=make_log_likelihood(n, counts_array),
                              #x0=np.array([(fix_c if main_allele == 'alt' else scaling[fix_c], 0.5)]),
                              x0=np.array([0.5]),
                              bounds=[(0, 1)])
    except ValueError:
        return 'NaN', 10
    w = x.x
    gof = calculate_gof(counts_array, w, r_array[fix_c])
    print(q_left, gof)
    return x, gof


def make_log_likelihood(n, counts_array):
    #print("target_made")

    def target(x):

        # print(r, p)
        # print("Counting likelihood")
        neg_bin_dens = make_negative_binom_density(r_array[fix_c], get_p(), x, right_most)
        return -1 * sum(counts_array[k] * np.log(neg_bin_dens[k])
                        for k in range(left_most, n) if counts_array[k] != 0)

    return target


def get_p():
    return 1 / (BAD + 1)


def calculate_gof(counts_array, w, r):
    observed = counts_array.copy()
    observed[:q_left] = 0
    norm = observed.sum()
    expected = make_negative_binom_density(r, get_p(), w, len(observed) - 1) * norm

    idxs = (observed != 0) & (expected != 0)
    df = idxs.sum() - 2

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
            precalc_params_path = parameters_path + 'NBweights_{}_BAD={:.1f}.npy'.format(fixed_allele, BAD)
            coefs_array = np.load(precalc_params_path)

            save_array = np.zeros((max(fix_c_array) + 1, 8), dtype=np.float_)
            save_array[:, :4] = coefs_array

            print(coefs_array.shape[0])

            max_idx = min(i for i in range(5, coefs_array.shape[0]) if coefs_array[i, 3] > 0.1)

            reg = LinearRegression().fit(X=np.array(range(5, max_idx)).reshape(-1, 1),
                                                            y=coefs_array[5: max_idx, 0],
                                                            sample_weight=1 / coefs_array[5: max_idx, 3])
            k = reg.coef_[0]
            b = reg.intercept_
            print(k, max_idx)

            r_array = [0] * 5 + [k * x for x in fix_c_array]

            filename = parameters_path + 'fixed_alt_bias_statistics_BAD={:.1f}.tsv'.format(BAD)
            stats = pd.read_table(filename)
            for allele in alleles:
                stats['{}_counts'.format(allele)] = stats['{}_counts'.format(allele)].astype(int)

            for fix_c in fix_c_array:
                stats_filtered = stats[stats['{}_counts'.format(fixed_allele)] == fix_c]
                counts, set_of_nonzero_n = make_counts(stats_filtered)
                if len(set_of_nonzero_n) == 0 or counts.sum() < max(set_of_nonzero_n) - 5:
                    continue
                print('made counts')
                print('Fix {}={}, k={}, max={}'.format(fixed_allele, fix_c, k, max_idx))
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
                    save_array[fix_c, 4] = k * fix_c + b
                    save_array[fix_c, 5] = weights.x
                    save_array[fix_c, 6] = weights.success
                    save_array[fix_c, 7] = gof

            np.save(os.path.expanduser(parameters_path + 'NBweights_step2_{}_BAD={:.1f}.npy'.format(
                fixed_allele, BAD)), save_array)
