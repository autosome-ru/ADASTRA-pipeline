import sys
import numpy as np
import pandas as pd
from math import ceil
from scipy import stats as st
from scripts.HELPERS.helpers import read_weights
from scripts.HELPERS.paths import get_ending

r_dict, w_dict, gof_dict = read_weights()


def count_p(ref_c, alt_c, BADs):
    N = len(ref_c)
    p_ref = np.zeros(N, dtype=np.float128)
    p_alt = np.zeros(N, dtype=np.float128)
    es_ref = np.zeros(N, dtype=np.float128)
    es_alt = np.zeros(N, dtype=np.float128)

    for i in range(N):
        if alt_c[i] >= ref_c[i] * BADs[i]:
            if ref_c[i] > 500:
                r = ref_c[i]
                r_ref = ref_c[i]
            else:
                r, gof = (r_dict['ref'][BADs[i]][ref_c[i]],
                          gof_dict['ref'][BADs[i]][ref_c[i]])
                r_ref = r_dict['alt'][BADs[i]][ref_c[i]]
                if r == 0:
                    r = ref_c[i]
                if r_ref == 0:
                    r_ref = ref_c[i]

            cdf = st.nbinom(r, 1 / (BADs[i] + 1)).cdf
            left_border = ceil(ref_c[i] * BADs[i] * r / r_ref)
            p_alt[i] = (1 - cdf(alt_c[i] - 1)) / (1 - cdf(left_border - 1))
            es_alt[i] = np.log2(alt_c[i] / left_border)
        else:
            p_alt[i] = np.nan
            es_alt[i] = np.nan

        #  _______________________________

        if ref_c[i] >= alt_c[i] * BADs[i]:
            if alt_c[i] > 500:
                r = alt_c[i]
            else:
                r, gof = (r_dict['alt'][BADs[i]][alt_c[i]],
                          gof_dict['alt'][BADs[i]][alt_c[i]])
                if r == 0:
                    r = alt_c[i]

            cdf = st.nbinom(r, 1 / (BADs[i] + 1)).cdf
            left_border = ceil(alt_c[i] * BADs[i])
            p_ref[i] = (1 - cdf(ref_c[i] - 1)) / (1 - cdf(left_border - 1))
            es_ref[i] = np.log2(ref_c[i] / left_border)
        else:
            p_ref[i] = np.nan
            es_ref[i] = np.nan
    return p_ref, p_alt, es_ref, es_alt


def test_pval():
    inp_str = ''
    while inp_str != 'q':
        inp_str = input('Enter ref_c alt_c BAD (enter q to exit):\n')
        try:
            refc, altc, BAD = [float(x) for x in inp_str.split()]
        except ValueError:
            continue
        p_ref, p_alt, es_ref, es_alt = count_p(np.array([refc], dtype=np.int_),
                                               np.array([altc], dtype=np.int_),
                                               np.array([BAD], dtype=np.float_))
        print('pval ref: {}\npval alt: {}\nes ref: {}\nes_alt: {}'.format(p_ref[0], p_alt[0], es_ref[0], es_alt[0]))


def main(base_path):
    table_BAD = base_path + get_ending("BAD")
    output = base_path + get_ending("p-value")
    print('Now counting P-value for {}'.format(table_BAD))
    df_with_BAD = pd.read_table(table_BAD)
    p_ref, p_alt, es_ref, es_alt = count_p(np.array(df_with_BAD["ref_read_counts"], dtype=np.int_),
                                           np.array(df_with_BAD["alt_read_counts"], dtype=np.int_),
                                           np.array(df_with_BAD["BAD"], dtype=np.float_))
    df_with_BAD['p_value_ref'] = p_ref
    df_with_BAD['p_value_alt'] = p_alt
    df_with_BAD['es_ref'] = es_ref
    df_with_BAD['es_alt'] = es_alt
    df_with_BAD.to_csv(output, sep="\t", index=False)


if __name__ == '__main__':
    main(sys.argv[1])
