import sys
import numpy as np
import pandas as pd

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import parameters_path
from scripts.HELPERS.helpers import states


def count_p(ref_c, alt_c, BAD):
    n = ref_c + alt_c
    w = precalc_data[BAD]['weights'][n]

    p_ref = (1 - w) * precalc_data[BAD]['binom_sum'][n, min(ref_c, alt_c)] +\
            w * precalc_data[BAD]['noise_sum_ref'][n, ref_c]
    p_alt = (1 - w) * precalc_data[BAD]['binom_sum'][n, max(ref_c, alt_c)] +\
            w * precalc_data[BAD]['noise_sum_alt'][n, ref_c]
    # assert 0 <= p_ref <= 1
    # assert 0 <= p_alt <= 1
    return p_ref, p_alt


if __name__ == '__main__':
    full_path = sys.argv[1]
    table_BAD = full_path + "_table_BADs.txt"
    output = full_path + "_table_p.txt"
    print('Now counting P-value for {}'.format(table_BAD))

    # loading precalculated density and noise weights
    precalc_data = {}
    for BAD in states:
        precalc_data[BAD] = {}
        filename = parameters_path + 'cover_bias_statistics_BAD={:.1f}.tsv'.format(BAD)
        precalc_data[BAD]['binom_sum'] = filename + '_binom_sum.precalc.npy'
        precalc_data[BAD]['noise_sum_ref'] = filename + '_noise_sum_ref.precalc.npy'
        precalc_data[BAD]['noise_sum_alt'] = filename + '_noise_sum_alt.precalc.npy'
        precalc_data[BAD]['weights'] = np.load(parameters_path + 'weights_BAD={:.1f}.npy'.format(BAD))

    df_with_BAD = pd.read_table(table_BAD)
    p_ref, p_alt = count_p(df_with_BAD["ref_read_counts"], df_with_BAD["alt_read_counts"],
                           df_with_BAD["BAD"])
    df_with_BAD['p_value_ref'] = p_ref
    df_with_BAD['p_value_alt'] = p_alt
    df_with_BAD.to_csv(output, sep="\t", index=False)
