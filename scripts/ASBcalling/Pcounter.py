import sys
import numpy as np
import pandas as pd
from scipy.stats import binom_test, binom

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import parameters_path
from scripts.HELPERS.helpers import states


def count_p(x, y, p):
    pv_ref = np.array(len(x), dtype=np.float_)
    pv_alt = np.array(len(x), dtype=np.float_)
    for alternative in ('alt', 'ref'):
        for i in range(len(x)):
            if alternative == 'ref':
                pv_ref[i] = (binom_test(x[i], x[i] + y[i], p, alternative) + binom_test(x[i], x[i] + y[i], 1 - p,
                                                                                        alternative)) / 2
            else:
                pv_alt[i] = (binom_test(x[i], x[i] + y[i], p, alternative) + binom_test(x[i], x[i] + y[i], 1 - p,
                                                                                        alternative)) / 2
    return pv_ref, pv_alt


# def count_p(ref_c, alt_c, BADs):
#     n = ref_c + alt_c
#     N = len(ref_c)
#     p_ref = np.zeros(N, dtype=np.float128)
#     p_alt = np.zeros(N, dtype=np.float128)
#
#     for BAD in np.unique(BADs):
#         # loading precalculated density and noise weights
#         precalc_data = {}
#         filename = parameters_path + 'cover_bias_statistics_BAD={:.1f}.tsv'.format(BAD)
#         precalc_data['binom_sum'] = np.load(filename + '_binom_sum.precalc.npy')
#         precalc_data['noise_sum_ref'] = np.load(filename + '_noise_sum_ref.precalc.npy')
#         precalc_data['noise_sum_alt'] = np.load(filename + '_noise_sum_alt.precalc.npy')
#         precalc_data['weights'] = np.load(parameters_path + 'weights_BAD={:.1f}.npy'.format(BAD))
#
#         idcs = np.where(BADs == BAD)
#         n_BAD = n[idcs]
#         ref_c_BAD = ref_c[idcs]
#         alt_c_BAD = alt_c[idcs]
#         w_BAD = precalc_data['weights'][n_BAD]
#
#         p_ref[idcs] = (1 - w_BAD) * precalc_data['binom_sum'][n_BAD, alt_c_BAD] + \
#                       w_BAD * precalc_data['noise_sum_ref'][n_BAD, ref_c_BAD]
#         p_alt[idcs] = (1 - w_BAD) * precalc_data['binom_sum'][n_BAD, ref_c_BAD] + \
#                       w_BAD * precalc_data['noise_sum_alt'][n_BAD, ref_c_BAD]
#     return p_ref, p_alt


if __name__ == '__main__':
    full_path = sys.argv[1]
    table_BAD = full_path + "_table_BADs.txt"
    output = full_path + "_table_p.txt"
    print('Now counting P-value for {}'.format(table_BAD))

    df_with_BAD = pd.read_table(table_BAD)
    p_ref, p_alt = count_p(np.array(df_with_BAD["ref_read_counts"], dtype=np.int_),
                           np.array(df_with_BAD["alt_read_counts"], dtype=np.int_),
                           np.array(df_with_BAD["BAD"], dtype=np.float_))
    df_with_BAD['p_value_ref'] = p_ref
    df_with_BAD['p_value_alt'] = p_alt
    df_with_BAD.to_csv(output, sep="\t", index=False)
