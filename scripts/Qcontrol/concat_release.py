import pandas as pd
import sys
import os
import json


def make_path(base_path, mode, is_dict=False):
    return os.path.join(base_path, mode + ('_DICTS' if is_dict else '_P-values'))


def read_dict(base_path, mode, tf):
    with open(os.path.join(make_path(base_path, mode, True), tf + '.json')) as js:
        d = json.loads(js.readline())
    return d


def get_id(row):
    return '\t'.join(map(str, [row['#chr'], row['pos'], row['ID'],
                               row['ref'], row['alt'],
                               '' if pd.isna(row['repeat_type']) else row['repeat_type']]))


def add_pvalue(row, d):
    snp_id = get_id(row)
    row['pval_ref'] = d[snp_id]["ref_pvalues"]
    row['pval_alt'] = d[snp_id]["alt_pvalues"]
    return row


def read_df(base_path, mode, tf, n_tr=50):
    df = pd.read_table(os.path.join(make_path(base_path, mode), tf + '.tsv'))
    df = df[df['n_aggregated'] >= n_tr]
    if not df.empty:
        df[mode] = tf
        d = read_dict(base_path, mode, tf)
        df = df.apply(lambda x: add_pvalue(x, d), axis=1)
    return df


def main(release_dir, modes=None):
    if modes is None:
        modes = ['CL']
    for mode in modes:
        dfs = []
        for file in os.listdir(make_path(release_dir, mode)):
            tf = os.path.splitext(file)[0]
            print('Now doing', tf)
            df = read_df(release_dir, mode, tf)
            if not df.empty:
                dfs.append(df)
        sum_df = pd.concat(dfs, ignore_index=True)
        print(sum_df)
        sum_df.to_csv('~/dan_conc_{}.tsv'.format(mode), sep='\t', index=False)


main(sys.argv[1], sys.argv[2].split(','))
