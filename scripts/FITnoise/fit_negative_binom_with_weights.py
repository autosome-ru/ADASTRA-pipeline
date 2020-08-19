import os
import numpy as np
import pandas as pd
from scipy import optimize
from scipy import stats as st
from scripts.HELPERS.paths_for_components import configs_path
from scripts.HELPERS.helpers import states


def make_negative_binom_density(r, p, w, size_of_counts, left_most, for_plot=False):
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


def make_scaled_counts(stats_pandas_dataframe, main_allele):
    try:
        max_cover_in_stats = max(stats_pandas_dataframe['{}_counts'.format(main_allele)])
    except ValueError:
        return [], set()
    counts_array = np.zeros(max_cover_in_stats + 1, dtype=np.int64)
    nonzero_set = set()

    for index, row in stats_pandas_dataframe.iterrows():
        k, SNP_counts = row['{}_counts'.format(main_allele)], row['counts']
        nonzero_set.add(k)
        counts_array[k] += SNP_counts
    return counts_array, nonzero_set


def fit_negative_binom(n, counts_array, fix_c, BAD, q_left, left_most):
    try:
        x = optimize.minimize(fun=make_log_likelihood(n, counts_array, BAD, left_most),
                              x0=np.array([fix_c, 0.5]),
                              bounds=[(0.00001, None), (0, 1)])
    except ValueError:
        return 'NaN', 10
    r, w = x.x
    return x, calculate_gof(counts_array, w, r, BAD, q_left, left_most)


def make_log_likelihood(n, counts_array, BAD, left_most):
    def target(x):
        r = x[0]
        w = x[1]
        neg_bin_dens = make_negative_binom_density(r, get_p(BAD), w, len(counts_array), left_most)
        return -1 * sum(counts_array[k] * np.log(neg_bin_dens[k])
                        for k in range(left_most, n) if counts_array[k] != 0)

    return target


def get_p(BAD):
    return 1 / (BAD + 1)


def calculate_gof(counts_array, w, r, BAD, q_left, left_most):
    # FIXME for BAD==1
    observed = counts_array.copy()
    observed[:q_left] = 0
    norm = observed.sum()
    expected = make_negative_binom_density(r, get_p(BAD), w, len(observed) - 1, left_most) * norm

    idxs = (observed != 0) & (expected != 0)
    df = idxs.sum() - 3

    stat = np.sum(observed[idxs] * np.log(observed[idxs] / expected[idxs])) * 2

    if norm <= 1:
        return None
    else:
        score = np.sqrt(max(stat - df, 0) / (df * (norm - 1)))

    return score


def main():
    for main_allele in ("alt", "ref"):
        fixed_allele = "ref" if main_allele == "alt" else "alt"
        alleles = ('ref', 'alt')

        fix_c_array = list(range(5, 501))

        for BAD in states:
            save_array = np.zeros((max(fix_c_array) + 1, 4), dtype=np.float_)
            filename = os.path.join(configs_path, 'bias_statistics_BAD={:.1f}.tsv'.format(BAD))
            stats = pd.read_table(filename)
            for allele in alleles:
                stats['{}_counts'.format(allele)] = stats['{}_counts'.format(allele)].astype(int)

            for fix_c in fix_c_array:
                stats_filtered = stats[stats['{}_counts'.format(fixed_allele)] == fix_c]
                counts, set_of_nonzero_n = make_scaled_counts(stats_filtered, main_allele)
                if len(set_of_nonzero_n) == 0 or counts.sum() < max(set_of_nonzero_n) - 5:
                    continue

                left_most = 5
                q_left = 5
                right_most = len(counts) - 1

                calculate_negative_binom = True
                if calculate_negative_binom:
                    weights, gof = fit_negative_binom(right_most, counts, fix_c, BAD, q_left, left_most)
                    save_array[fix_c, :2] = weights.x
                    save_array[fix_c, 2] = weights.success
                    save_array[fix_c, 3] = gof

            np.save(os.path.join(configs_path, 'NBweights_{}_BAD={:.1f}'.format(fixed_allele, BAD)), save_array)


if __name__ == '__main__':
    main()
