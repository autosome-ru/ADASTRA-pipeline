import pandas as pd
import os
from scripts.HELPERS.paths_for_components import master_list_path, parallel_parameters_path

from scripts.HELPERS.helpers import dtype_dict


def main():
    master_df = pd.read_table(master_list_path, dtype=dtype_dict)
    if 'DOWNLOAD_PATH' not in master_df.columns:
        master_df['DOWNLOAD_PATH'] = [''] * (len(master_df.index))
    else:
        master_df = master_df[(master_df['DOWNLOAD_PATH'] != 'None') & (master_df['DOWNLOAD_PATH'].notna())]
    master_df = master_df.sort_values(
        by=['ALIGNED_READS'],
        axis=0,
        ascending=False,
        na_position='first')
    master_df = master_df[['#EXP', 'ALIGNS', 'READS', 'DOWNLOAD_PATH']]
    master_df.to_csv(os.path.join(parallel_parameters_path, 'sorted_master_list.tsv'), sep='\t', index=False, header=False)


if __name__ == '__main__':
    main()
