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
        if BADs[i] != 1:
            p_ref[i] = 'NaN'
            p_alt[i] = 'NaN'
        else:
            fixed = min(ref_c[i], alt_c[i])
            main = max(ref_c[i], alt_c[i])
            cdf = (st.nbinom(fixed, 1 / (BADs[i] + 1))).cdf
            pv = (1 - cdf(main)) / (1 - cdf(4))
            if ref_c[i] > alt_c[i]:
                p_alt[i] = 'NaN'
                p_ref[i] = pv
            elif ref_c[i] < alt_c[i]:
                p_alt[i] = pv
                p_ref[i] = 'NaN'
            else:
                p_alt[i] = 'NaN'
                p_ref[i] = 'NaN'

    return p_ref, p_alt


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
