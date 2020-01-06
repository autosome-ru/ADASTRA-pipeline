import sys
import numpy as np
import pandas as pd
from scipy import stats as st

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import parameters_path
from scripts.HELPERS.helpers import states


def count_p(ref_c, alt_c, BADs):
    N = len(ref_c)
    p_ref = np.zeros(N, dtype=np.float128)
    p_alt = np.zeros(N, dtype=np.float128)

    for i in range(N):

        if ref_c[i] > 500:
            r = ref_c[i]
            w = 0.5
        else:
            r, w = weights[BADs[i]]['ref'][ref_c[i], :2]
            if r == 0:
                r = ref_c[i]
                w = 0.5
        cdf1 = (st.nbinom(r, 1 / (BADs[i] + 1))).cdf
        cdf2 = (st.nbinom(r, BADs[i] / (BADs[i] + 1))).cdf
        cdf = lambda x: w * cdf1(x) + (1 - w) * cdf2(x)
        p_alt[i] = (1 - cdf(alt_c[i])) / (1 - cdf(4))

        if alt_c[i] > 500:
            r = alt_c[i]
            w = 0.5
        else:
            r, w = weights[BADs[i]]['alt'][alt_c[i], :2]
            if r == 0:
                r = alt_c[i]
                w = 0.5
        cdf1 = (st.nbinom(r, 1 / (BADs[i] + 1))).cdf
        cdf2 = (st.nbinom(r, BADs[i] / (BADs[i] + 1))).cdf
        cdf = lambda x: w * cdf1(x) + (1 - w) * cdf2(x)
        p_ref[i] = (1 - cdf(ref_c[i])) / (1 - cdf(4))
    return p_ref, p_alt


if __name__ == '__main__':
    full_path = sys.argv[1]
    table_BAD = full_path + "_table_BADs.txt"
    output = full_path + "_table_p.txt"
    print('Now counting P-value for {}'.format(table_BAD))

    weights = {}
    for BAD in states:
        weights[BAD] = dict()
        weights[BAD]['ref'] = np.load(parameters_path + 'NBweights_ref_BAD={:.1f}'.format(BAD))
        weights[BAD]['alt'] = np.load(parameters_path + 'NBweights_alt_BAD={:.1f}'.format(BAD))

    df_with_BAD = pd.read_table(table_BAD)
    p_ref, p_alt = count_p(np.array(df_with_BAD["ref_read_counts"], dtype=np.int_),
                           np.array(df_with_BAD["alt_read_counts"], dtype=np.int_),
                           np.array(df_with_BAD["BAD"], dtype=np.float_))
    df_with_BAD['p_value_ref'] = p_ref
    df_with_BAD['p_value_alt'] = p_alt
    df_with_BAD.to_csv(output, sep="\t", index=False)
