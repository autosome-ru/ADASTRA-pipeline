import os
import pandas as pd
from scripts.HELPERS.paths import get_result_dir_path, get_result_table_path
from scripts.SARUSannotation.extract_sarus_data import main as main_extract


def main():
    all_asbs_df = None
    for obj in 'TF', 'CL':
        for tf_file in os.listdir(get_result_dir_path(obj)):
            tf_name = os.path.splitext(tf_file)[0]
            tf_df = pd.read_table(get_result_table_path(obj, tf_name))
            tf_df = tf_df[tf_df.columns.drop(list(tf_df.filter(regex='motif')))]
            tf_df = tf_df[(tf_df['fdrp_bh_ref'] <= 0.25) | (tf_df['fdrp_bh_alt'] <= 0.25)]
            if tf_df.empty:
                continue
            tf_df['unique'] = '{}@{}{}'.format(tf_df['#chr'],
                                               tf_df['pos'],
                                               tf_df['alt'])
            if all_asbs_df is None:
                all_asbs_df = tf_df
            else:
                all_asbs_df = all_asbs_df.append(tf_df, ignore_index=True)
    main_extract('all_tfs', 25, all_asbs_df.sort_values(by=['pos', '#chr']))
