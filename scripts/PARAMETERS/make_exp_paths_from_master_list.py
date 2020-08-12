import json
import pandas as pd
import os.path
import sys

from scripts.HELPERS.paths import get_ending, create_path_from_master_list_df
from scripts.HELPERS.paths_for_components import parallel_parameters_path, badmaps_dict_path, master_list_path

out_path = os.path.join(parallel_parameters_path, 'exp_paths.cfg')


if __name__ == "__main__":
    if len(sys.argv) > 1:
        for_what = sys.argv[1]
    else:
        for_what = 'badmaps'
    master_df = pd.read_table(master_list_path)
    master_df = master_df[master_df['EXP_TYPE'] != 'chip_control']
    master_df['path'] = master_df.apply(create_path_from_master_list_df)
    master_df = master_df[master_df['path'].apply(
        lambda x: os.path.isfile(x + get_ending('vcf')))]
    if for_what == 'badmaps':
        master_df = master_df[master_df['path'].apply(
            lambda x: os.path.isfile(x + get_ending('annotation')))]
        master_df['path'].to_csv(out_path, sep='\t', index=False, header=False)
    elif for_what == 'annotation':
        master_df[['path', 'PEAKS']].to_csv(out_path, sep='\t', index=False, header=False)
