import pandas as pd
import os.path
import sys

from scripts.HELPERS.helpers import dtype_dict
from scripts.HELPERS.paths import get_ending, create_path_from_master_list_df
from scripts.HELPERS.paths_for_components import parallel_parameters_path, badmaps_dict_path, master_list_path

out_path = os.path.join(parallel_parameters_path, 'exp_paths.cfg')


def main(for_what):
    master_df = pd.read_table(master_list_path, dtype=dtype_dict)
    master_df = master_df[~master_df['EXP_TYPE'].isin(['chip_control', 'chipexo_control'])]
    master_df['path'] = master_df.apply(create_path_from_master_list_df, axis=1)
    master_df = master_df[master_df['path'].apply(
        lambda x: os.path.isfile(x + get_ending('vcf')))]
    if for_what == 'badmaps':
        master_df = master_df[master_df['path'].apply(
            lambda x: os.path.isfile(x + get_ending('annotation')))]
        master_df['path'].to_csv(out_path, sep='\t', index=False, header=False)
    elif for_what == 'annotation':
        master_df[['path', 'PEAKS']].to_csv(out_path, sep='\t', index=False, header=False)


if __name__ == '__main__':
    main(sys.argv[1])
