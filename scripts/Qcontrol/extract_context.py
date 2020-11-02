import os
import pandas as pd
from scripts.HELPERS.paths import get_result_dir_path, get_result_table_path
from scripts.SARUSannotation.extract_sarus_data import main as main_extract


def main():
    all_asbs_df = None
    for tf_file in os.listdir(get_result_dir_path('TF')):
        tf_name = os.path.splitext(tf_file)[0]
        tf_df = pd.read_table(get_result_table_path('TF', tf_name))
        tf_df = tf_df[(tf_df['fdrp_bh_ref'] <= 0.05) | (tf_df['fdrp_bh_alt'] <= 0.05)]
        tf_df['unique'] = '{}@{}{}'.format(all_asbs_df['#chr'],
                                           all_asbs_df['pos'],
                                           all_asbs_df['alt'])
        if all_asbs_df is None:
            all_asbs_df = tf_df
        else:
            all_asbs_df = all_asbs_df.append(tf_df, ignore_index=True)
        all_asbs_df.drop_duplicates(subset='unique')
        print(len(all_asbs_df.index))
    main_extract('all_tfs', 25, all_asbs_df)
