import json

import pandas as pd
import os.path
import sys

from scripts.HELPERS.helpers import dtype_dict, make_reverse_dict, is_valid, get_merged_badmaps_dict_path
from scripts.HELPERS.paths import get_ending, create_path_from_master_list_df
from scripts.HELPERS.paths_for_components import parallel_parameters_path, master_list_path

out_path = os.path.join(parallel_parameters_path, 'exp_paths.cfg')


def main(for_what, remade=True):
    master_df = pd.read_table(master_list_path, dtype=dtype_dict)
    master_df = master_df[~master_df['EXP_TYPE'].isin(['chip_control', 'chipexo_control'])]
    master_df['path'] = master_df.apply(create_path_from_master_list_df, axis=1)
    if for_what == 'badmaps':
        with open(get_merged_badmaps_dict_path(remade=remade), "r") as read_file:
            d = json.load(read_file)
            rev_d = make_reverse_dict(d)
        df_filter = master_df.apply(
            lambda row: os.path.isfile(create_path_from_master_list_df(row, 'BAD'))
                        and is_valid(row['path'], rev_d, remade=remade), axis=1)
        print(len(df_filter.sum()))
        master_df = master_df[df_filter]
        master_df['path'].to_csv(out_path, sep='\t', index=False, header=False)
    elif for_what == 'annotation':
        master_df[['path', 'PEAKS']].to_csv(out_path, sep='\t', index=False, header=False)


if __name__ == '__main__':
    main(sys.argv[1])
