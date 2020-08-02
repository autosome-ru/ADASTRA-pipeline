import sys
import numpy as np
import pandas as pd
from scipy import stats as st
from scripts.HELPERS.helpers import read_weights


def count_p(ref_c, alt_c, BADs):
    N = len(ref_c)
    p_ref = np.zeros(N, dtype=np.float128)
    p_alt = np.zeros(N, dtype=np.float128)
    es_ref = np.zeros(N, dtype=np.float128)
    es_alt = np.zeros(N, dtype=np.float128)

    for i in range(N):

        if ref_c[i] > 500:
            r = ref_c[i]
            w = 1
        else:
            r, w, gof = (r_dict['ref'][BADs[i]][ref_c[i]],
                         w_dict['ref'][BADs[i]][ref_c[i]],
                         gof_dict['ref'][BADs[i]][ref_c[i]])
            if r == 0:
                r = ref_c[i]
                w = 1
        dist1 = st.nbinom(r, 1 / (BADs[i] + 1))
        dist2 = st.nbinom(r, BADs[i] / (BADs[i] + 1))
        cdf1 = dist1.cdf
        cdf2 = dist2.cdf
        pmf1 = dist1.pmf
        pmf2 = dist2.pmf
        pmf = lambda x: w * pmf1(x) + (1 - w) * pmf2(x)
        cdf = lambda x: w * cdf1(x) + (1 - w) * cdf2(x)
        p_alt[i] = (1 - cdf(alt_c[i] - 1)) / (1 - cdf(4))
        if p_alt[i] != 1:
            E_alt = (r * (BADs[i] * w + (1 - w) / BADs[i]) - sum(i * pmf(i) for i in range(5))) / (
                    1 - cdf(4))
            es_alt[i] = np.log(alt_c[i] / E_alt)
        else:
            es_alt[i] = 'NaN'

        if alt_c[i] > 500:
            r = alt_c[i]
            w = 1
        else:
            r, w, gof = (r_dict['alt'][BADs[i]][alt_c[i]],
                         w_dict['alt'][BADs[i]][alt_c[i]],
                         gof_dict['alt'][BADs[i]][alt_c[i]])
            if r == 0:
                r = alt_c[i]
                w = 1
        dist1 = st.nbinom(r, 1 / (BADs[i] + 1))
        dist2 = st.nbinom(r, BADs[i] / (BADs[i] + 1))
        cdf1 = dist1.cdf
        cdf2 = dist2.cdf
        pmf1 = dist1.pmf
        pmf2 = dist2.pmf
        pmf = lambda x: w * pmf1(x) + (1 - w) * pmf2(x)
        cdf = lambda x: w * cdf1(x) + (1 - w) * cdf2(x)
        p_ref[i] = (1 - cdf(ref_c[i] - 1)) / (1 - cdf(4))
        if p_ref[i] != 1:
            E_ref = (r * (BADs[i] * w + (1 - w) / BADs[i]) - sum(i * pmf(i) for i in range(5))) / (
                    1 - cdf(4))
            es_ref[i] = np.log(ref_c[i] / E_ref)
        else:
            es_ref[i] = 'NaN'
    return p_ref, p_alt, es_ref, es_alt


if __name__ == '__main__':
    full_path = sys.argv[1]
    table_BAD = full_path + "_table_BADs.txt"
    output = full_path + "_table_p.txt"
    print('Now counting P-value for {}'.format(table_BAD))

    r_dict, w_dict, gof_dict = read_weights()

    df_with_BAD = pd.read_table(table_BAD)
    p_ref, p_alt, es_ref, es_alt = count_p(np.array(df_with_BAD["ref_read_counts"], dtype=np.int_),
                           np.array(df_with_BAD["alt_read_counts"], dtype=np.int_),
                           np.array(df_with_BAD["BAD"], dtype=np.float_))
    df_with_BAD['p_value_ref'] = p_ref
    df_with_BAD['p_value_alt'] = p_alt
    df_with_BAD['es_ref'] = es_ref
    df_with_BAD['es_alt'] = es_alt
    df_with_BAD.to_csv(output, sep="\t", index=False)
