import os
import pandas as pd
import sys
import numpy as np
from tqdm import tqdm

results_path = sys.argv[1]


def get_result_dir_path(what_for):
    return os.path.join(results_path, '{}_P-values'.format(what_for))


def calc_conc(row):
    row['motif_conc'] = get_concordance(row['fdrp_bh_ref'], row['fdrp_bh_alt'],
                           row['motif_fc'], row['motif_log_pref'],
                           row['motif_log_palt'])
    return row


def get_concordance(p_val_ref, p_val_alt, motif_fc, motif_pval_ref, motif_pval_alt):
    if not pd.isna(p_val_ref) and not pd.isna(p_val_alt):
        log_pv = np.log10(min(p_val_ref, p_val_alt)) * np.sign(p_val_alt - p_val_ref)
        if abs(log_pv) >= -np.log10(0.25):
            if max(motif_pval_ref, motif_pval_alt) >= -np.log10(0.0005) and motif_fc != 0:
                result = "Weak " if abs(motif_fc) < 2 else ""
                if motif_fc * log_pv > 0:
                    result += 'Concordant'
                elif motif_fc * log_pv < 0:
                    result += 'Discordant'
                return result
            else:
                return "No Hit"
    return None


def get_result_table_path(what_for, string):
    return os.path.join(get_result_dir_path(what_for), '{}.tsv'.format(string))


def main():
    for obj in 'TF', 'CL':
        tf_files = os.listdir(get_result_dir_path(obj))
        print(f'Doing {obj}')
        for tf_file in tqdm(tf_files):
            tf_name = os.path.splitext(tf_file)[0]
            tf_df = pd.read_table(get_result_table_path(obj, tf_name),
                                  dtype=None if obj == 'CL' else {'motif_conc': str})
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
            if obj == 'TF':
                tf_df = tf_df.apply(calc_conc, axis=1)
            tf_df = tf_df[(~tf_df['fdrp_bh_ref'].isna()) & (~tf_df['fdrp_bh_alt'].isna())]
            if tf_df.empty:
                continue
            tf_df.to_csv(os.path.join('~/release_dump/', obj, tf_file),
                         sep='\t', index=False)

main()
