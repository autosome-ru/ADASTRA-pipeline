import sys
import os.path
import json
import pandas as pd
import numpy as np


sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths_for_components import parameters_path, cl_dict_path
from scripts.HELPERS.helpers import remove_punctuation


def collectFixedAltStatistics(key_name=None, BAD=None, suffix=''):
    with open(cl_dict_path, "r") as read_file:
        cell_lines_dict = json.loads(read_file.readline())
    out_t = None

    if key_name:
        key_name = {remove_punctuation(key) for key in key_name}

    for key in cell_lines_dict:
        if key_name is not None:
            if key not in key_name:  # <------
                continue
        for align_path in cell_lines_dict[key]:
            bad_table_path = align_path.replace("_table_p", "_table_BADs")
            if not os.path.isfile(bad_table_path):
                continue
            print(bad_table_path)
            df = pd.read_table(bad_table_path)
            if df.empty:
                continue
            df = df[df['ID'].startswith('rs')]
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
    with open(parameters_path + 'fixed_alt_bias_statistics_BAD={:.1f}{}.tsv'.format(BAD if BAD else 0, suffix), 'w') as out:
        out_t.to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    esc = {'H1 (embryonic stem cells)',
             'BGO3 (embryonic stem cells)',
             'HUES64 (embryonic stem cells)',
             'CyT49 (embryonic stem cells)',
             'H9 (embryonic stem cells)',
             'VAL-3 (embryonic stem cells)',
             'WA09 (embryonic stem cells)',
             'embryonic stem cells',
             'human embryonic stem cells, H1 (WA01)'}
    hct116 = {"HCT-116 (colon carcinoma)"}
    k562 = {'K562 (myelogenous leukemia)'}
    for BAD in [None, 1, 2, 3, 4, 5, 6, 4/3, 5/2, 3/2]:
        collectFixedAltStatistics(BAD=BAD, key_name=None, suffix='')
        collectFixedAltStatistics(BAD=BAD, key_name=esc, suffix='_esc')
        collectFixedAltStatistics(BAD=BAD, key_name=hct116, suffix='_hct116')
        collectFixedAltStatistics(BAD=BAD, key_name=k562, suffix='_k562')
