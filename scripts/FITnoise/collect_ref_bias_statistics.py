import json
import os.path
import pandas as pd
from scripts.HELPERS.paths_for_components import master_list_path
from scripts.HELPERS.helpers import make_reverse_dict, is_valid, segmentation_states, dtype_dict, \
    get_merged_badmaps_dict_path
from scripts.HELPERS.paths import create_path_from_master_list_df, \
    create_neg_bin_stats_path_function


def collect_fixed_alt_statistics(master_df, key_name=None, BAD=None, suffix='', remade=True):
    out_t = None
    with open(get_merged_badmaps_dict_path(remade=remade), "r") as read_file:
        d = json.load(read_file)
        rev_d = make_reverse_dict(d)
    for index, row in master_df.iterrows():
        if key_name is not None:
            if row['CELLS'] not in key_name:  # <------
                continue
        base_path = create_path_from_master_list_df(row)
        if not is_valid(base_path, rev_d, remade=remade):
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
            out_t = out_t.groupby(['ref_counts', 'alt_counts']).size().reset_index(name='counts')
            out_t.fillna(0, inplace=True)
        else:
            tmp_df = pd.DataFrame()
            tmp_df['alt_counts'] = sum_df['alt_read_counts']
            tmp_df['ref_counts'] = sum_df['ref_read_counts']
            tmp_df = tmp_df.groupby(['ref_counts', 'alt_counts']).size().reset_index(name='counts')
            tmp_df.fillna(0, inplace=True)
            out_t = out_t.append(tmp_df).groupby(['ref_counts', 'alt_counts'], as_index=False).sum()
    if out_t is None:
        return
    out_t.to_csv(create_neg_bin_stats_path_function(BAD, suffix), sep="\t", index=False)


def main(cell_line=None, suffix='', in_stats=False, remade=True):
    master_df = pd.read_table(master_list_path, dtype=dtype_dict)
    master_df = master_df[~master_df['EXP_TYPE'].isin(['chip_control', 'chipexo_control'])]
    if in_stats:
        if not cell_line:
            collect_fixed_alt_statistics(master_df, BAD=None, remade=remade)
        else:
            for bad in [None] + segmentation_states:
                collect_fixed_alt_statistics(master_df, BAD=bad, key_name=cell_line, suffix=suffix, remade=remade)
    else:
        for bad in segmentation_states:
            collect_fixed_alt_statistics(master_df, BAD=bad, remade=remade)

