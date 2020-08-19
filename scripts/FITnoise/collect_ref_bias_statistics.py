import os.path
import pandas as pd

from scripts.HELPERS.helpers import dtype_dict
from scripts.HELPERS.paths import create_path_from_master_list_df
from scripts.HELPERS.paths_for_components import configs_path, master_list_path


def collect_fixed_alt_statistics(key_name=None, BAD=None, suffix=''):
    out_t = None
    master_df = pd.read_table(master_list_path, dtype=dtype_dict)
    master_df = master_df[master_df['EXP_TYPE'] != 'chip_control']
    for index, row in master_df.iterrows():
        if key_name is not None:
            if row['CELLS'] not in key_name:  # <------
                continue
        bad_table_path = create_path_from_master_list_df(row, 'BAD')
        if not os.path.isfile(bad_table_path):
            continue
        df = pd.read_table(bad_table_path)
        if df.empty:
            continue
        if BAD is not None:
            sum_df = df[df['BAD'] == BAD][['ref_read_counts', 'alt_read_counts']]  # <------
        else:
            sum_df = df[['ref_read_counts', 'alt_read_counts']]

        if out_t is None:
            out_t = pd.DataFrame()
            out_t['alt_counts'] = sum_df['alt_read_counts']
            out_t['ref_counts'] = sum_df['ref_read_counts']
            out_t = out_t.groupby(['alt_counts', 'ref_counts']).size().reset_index(name='counts')
            out_t.fillna(0, inplace=True)
        else:
            tmp_df = pd.DataFrame()
            tmp_df['alt_counts'] = sum_df['alt_read_counts']
            tmp_df['ref_counts'] = sum_df['ref_read_counts']
            tmp_df = tmp_df.groupby(['alt_counts', 'ref_counts']).size().reset_index(name='counts')
            tmp_df.fillna(0, inplace=True)
            out_t = out_t.append(tmp_df).groupby(['alt_counts', 'ref_counts'], as_index=False).sum()
    if out_t is None:
        return
    with open(os.path.join(configs_path, 'bias_statistics_BAD={:.1f}{}.tsv'.format(
            BAD if BAD else 0, suffix)), 'w') as out:
        out_t.to_csv(out, sep="\t", index=False)


def main():
    for bad in [1, 2, 3, 4, 5, 6, 4/3, 5/2, 3/2]:
        collect_fixed_alt_statistics(BAD=bad, key_name=None, suffix='')


if __name__ == "__main__":
    for BAD in [None, 1, 2, 3, 4, 5, 6, 4/3, 5/2, 3/2]:
        collect_fixed_alt_statistics(BAD=BAD, key_name=None, suffix='')
