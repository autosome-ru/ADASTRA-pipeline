import sys
import pandas as pd

sys.path.insert(1, "/home/abramov/ASB-Project")


def count_p(x, n, p):
    return None, None


if __name__ == '__main__':
    full_path = sys.argv[1]
    table_BAD = full_path + "_table_BADs.txt"
    output = full_path + "_table_p.txt"
    print('Now counting P-value for {}'.format(table_BAD))

    df_with_BAD = pd.read_table(table_BAD)
    p_ref, p_alt = count_p(df_with_BAD["ref_read_counts"], df_with_BAD["alt_read_counts"],
                           df_with_BAD["BAD"])
    df_with_BAD['p_value_ref'] = p_ref
    df_with_BAD['p_value_alt'] = p_alt
    df_with_BAD.to_csv(output, sep="\t", index=False)
