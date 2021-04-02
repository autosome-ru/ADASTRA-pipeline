import sys
import numpy as np
import pandas as pd
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
        if (i+1) % 1000 == 0:
            print(i+1)
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
        w = w * pmf1(alt_c[i]) / (w * pmf1(alt_c[i]) + (1-w) * pmf2(alt_c[i]))  # posterior w

        cdf = lambda x: w * cdf1(x) + (1 - w) * cdf2(x)
        pmf = lambda x: w * pmf1(x) + (1 - w) * pmf2(x)
        p_alt[i] = (1 - cdf(alt_c[i] - 1)) / (1 - cdf(4))
        if p_alt[i] != 1:
            E_alt = (r * (BADs[i] * w + (1 - w) / BADs[i]) - sum(i * pmf(i) for i in range(5))) / (
                    1 - cdf(4))
            es_alt[i] = np.log(alt_c[i] / E_alt)
        else:
            es_alt[i] = 'NaN'

        # ----------------------

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

        w = w * pmf1(ref_c[i]) / (w * pmf1(ref_c[i]) + (1-w) * pmf2(ref_c[i]))  # posterior w

        cdf = lambda x: w * cdf1(x) + (1 - w) * cdf2(x)
        pmf = lambda x: w * pmf1(x) + (1 - w) * pmf2(x)
        p_ref[i] = (1 - cdf(ref_c[i] - 1)) / (1 - cdf(4))
        if p_ref[i] != 1:
            E_ref = (r * (BADs[i] * w + (1 - w) / BADs[i]) - sum(i * pmf(i) for i in range(5))) / (
                    1 - cdf(4))
            es_ref[i] = np.log(ref_c[i] / E_ref)
        else:
            es_ref[i] = 'NaN'
    return p_ref, p_alt, es_ref, es_alt


def count_p_adjusted(ref_c, alt_c, BADs):
    N = len(ref_c)
    p_ref = np.zeros(N, dtype=np.float128)
    p_alt = np.zeros(N, dtype=np.float128)
    p_ref_bayes = np.zeros(N, dtype=np.float128)
    p_alt_bayes = np.zeros(N, dtype=np.float128)
    p_ref_likelihood = np.zeros(N, dtype=np.float128)
    p_alt_likelihood = np.zeros(N, dtype=np.float128)

    for i in range(N):
        if (i+1) % 1000 == 0:
            print(i+1)
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
        w_bayes = w * pmf1(alt_c[i]) / (w * pmf1(alt_c[i]) + (1-w) * pmf2(alt_c[i]))
        w_likelihood = (1 - w) * pmf1(alt_c[i]) / pmf2(alt_c[i])

        cdf = lambda x: w * cdf1(x) + (1 - w) * cdf2(x)
        p_alt[i] = (1 - cdf(alt_c[i] - 1)) / (1 - cdf(4))

        cdf = lambda x: w_bayes * cdf1(x) + (1 - w_bayes) * cdf2(x)
        p_alt_bayes[i] = (1 - cdf(alt_c[i] - 1)) / (1 - cdf(4))

        cdf = lambda x: w_likelihood * cdf1(x) + (1 - w_likelihood) * cdf2(x)
        p_alt_likelihood[i] = (1 - cdf(alt_c[i] - 1)) / (1 - cdf(4))

        # ----------------------

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
        w_bayes = w * pmf1(ref_c[i]) / (w * pmf1(ref_c[i]) + (1-w) * pmf2(ref_c[i]))
        w_likelihood = (1 - w) * pmf1(ref_c[i]) / pmf2(ref_c[i])

        cdf = lambda x: w * cdf1(x) + (1 - w) * cdf2(x)
        p_ref[i] = (1 - cdf(ref_c[i] - 1)) / (1 - cdf(4))

        cdf = lambda x: w_bayes * cdf1(x) + (1 - w_bayes) * cdf2(x)
        p_ref_bayes[i] = (1 - cdf(ref_c[i] - 1)) / (1 - cdf(4))

        cdf = lambda x: w_likelihood * cdf1(x) + (1 - w_likelihood) * cdf2(x)
        p_ref_likelihood[i] = (1 - cdf(ref_c[i] - 1)) / (1 - cdf(4))

    return (p_ref, p_alt,
           p_ref_bayes, p_alt_bayes,
           p_ref_likelihood, p_alt_likelihood)


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


def manual(exp, aligns):
    table_BAD = '/home/abramov/AlignmentsChip/{}/{}'.format(exp, aligns) + get_ending("BAD")
    output = '/home/abramov/test_K562_weighted_p/{}_{}'.format(exp, aligns) + get_ending("p-value")
    print('Now counting P-value for {}'.format(table_BAD))
    df_with_BAD = pd.read_table(table_BAD)
    # df_with_BAD = df_with_BAD[df_with_BAD['#chr'] == 'chr2']
    print(len(df_with_BAD.index))
    (p_ref, p_alt,
     p_ref_bayes, p_alt_bayes,
     p_ref_likelihood, p_alt_likelihood) = count_p_adjusted(np.array(df_with_BAD["ref_read_counts"], dtype=np.int_),
                                           np.array(df_with_BAD["alt_read_counts"], dtype=np.int_),
                                           np.array(df_with_BAD["BAD"], dtype=np.float_))
    df_with_BAD['p_value_ref'] = p_ref
    df_with_BAD['p_value_alt'] = p_alt
    df_with_BAD['p_value_ref_bayes'] = p_ref_bayes
    df_with_BAD['p_value_alt_bayes'] = p_alt_bayes
    df_with_BAD['p_value_ref_likelihood'] = p_ref_likelihood
    df_with_BAD['p_value_alt_likelihood'] = p_alt_likelihood

    print('i dump..')
    df_with_BAD.to_csv(output, sep="\t", index=False)
    print('i dump!')


if __name__ == '__main__':
    manual(sys.argv[1], sys.argv[2])
