import os
import pandas as pd
import sys

results_path = sys.argv[1]


def get_result_dir_path(what_for):
    return os.path.join(results_path, '{}_P-values'.format(what_for))


def get_result_table_path(what_for, string):
    return os.path.join(get_result_dir_path(what_for), '{}.tsv'.format(string))


def main():
    for obj in 'TF', 'CL':
        for tf_file in os.listdir(get_result_dir_path(obj)):
            tf_name = os.path.splitext(tf_file)[0]
            tf_df = pd.read_table(get_result_table_path(obj, tf_name))
            tf_df = tf_df[tf_df.columns.drop(['logitp_ref',
                                              'logitp_alt',
                                              'median_cover',
                                              'max_cover',
                                              'min_cover',
                                              'refc_mostsig_ref',
                                              'altc_mostsig_ref',
                                              'BAD_mostsig_ref',
                                              'es_mostsig_ref',
                                              'p_mostsig_ref',
                                              'refc_mostsig_alt',
                                              'altc_mostsig_alt',
                                              'BAD_mostsig_alt',
                                              'es_mostsig_alt',
                                              'p_mostsig_alt'
                                              ])]
            tf_df = tf_df[(~tf_df['fdrp_bh_ref'].isna()) & (~tf_df['fdrp_bh_alt'].isna())]
            if tf_df.empty:
                continue
            tf_df.to_csv(os.path.join('~/release_dump/', obj, tf_file),
                         sep='\t', index=False)

main()
