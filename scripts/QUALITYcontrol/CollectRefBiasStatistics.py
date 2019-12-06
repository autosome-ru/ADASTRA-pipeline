import sys
import os.path
import json
import pandas as pd

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import cl_dict_path, parameters_path

def collectRefAltStatistics(key_name=None, BAD=None):
    with open(cl_dict_path, "r") as read_file:
        cell_lines_dict = json.loads(read_file.readline())
    out_t = None

    for key in cell_lines_dict:
        if key_name is not None:
            if key not in key_name:  # <------
                continue
        for align_path in cell_lines_dict[key]:
            if not os.path.isfile(align_path):
                continue
            print(align_path)
            df = pd.read_table(align_path)
            if df.empty:
                continue
            if BAD is not None:
                sum_df = df[df['BAD'] == BAD][['ref_read_counts', 'alt_read_counts']] # <------
            else:
                sum_df = df[['ref_read_counts', 'alt_read_counts']]

            if out_t is None:
                out_t = pd.DataFrame()
                out_t['ref'] = sum_df['ref_read_counts'].value_counts()
                out_t['alt'] = sum_df['alt_read_counts'].value_counts()
                out_t['allele_reads'] = out_t.index
                out_t = out_t.reset_index(drop=True)
                out_t.fillna(0, inplace=True)
            else:
                tmp_df = pd.DataFrame()
                tmp_df['ref'] = sum_df['ref_read_counts'].value_counts()
                tmp_df['alt'] = sum_df['alt_read_counts'].value_counts()
                tmp_df['allele_reads'] = tmp_df.index
                tmp_df = tmp_df.reset_index(drop=True)
                tmp_df.fillna(0, inplace=True)
                out_t = out_t.append(tmp_df).groupby('allele_reads', as_index=False).sum()
    if out_t is None:
        return
    with open(parameters_path + 'bias_statistics.tsv', 'w') as out:
        out_t.to_csv(out, sep="\t")


def collectCoverStatistics(key_name=None, BAD=None):
    with open(cl_dict_path, "r") as read_file:
        cell_lines_dict = json.loads(read_file.readline())
    out_t = None

    for key in cell_lines_dict:
        if key_name is not None:
            if key not in key_name:  # <------
                continue
        for align_path in cell_lines_dict[key]:
            if not os.path.isfile(align_path):
                continue
            print(align_path)
            df = pd.read_table(align_path)
            if df.empty:
                continue
            if BAD is not None:
                sum_df = df[df['BAD'] >= BAD][['ref_read_counts', 'alt_read_counts']]  # <------
            else:
                sum_df = df[['ref_read_counts', 'alt_read_counts']]

            if out_t is None:
                out_t = pd.DataFrame()
                out_t['cover'] = (sum_df['ref_read_counts'] + sum_df['alt_read_counts'])
                out_t['ref_counts'] = sum_df['ref_read_counts']
                out_t = out_t.groupby(['cover', 'ref_counts']).size().reset_index(name='counts')
                out_t.fillna(0, inplace=True)
            else:
                tmp_df = pd.DataFrame()
                tmp_df['cover'] = (sum_df['ref_read_counts'] + sum_df['alt_read_counts'])
                tmp_df['ref_counts'] = sum_df['ref_read_counts']
                tmp_df = tmp_df.groupby(['cover', 'ref_counts']).size().reset_index(name='counts')
                tmp_df.fillna(0, inplace=True)
                out_t = out_t.append(tmp_df).groupby(['cover', 'ref_counts'], as_index=False).sum()
    if out_t is None:
        return
    with open(parameters_path + 'cover_bias_statistics.tsv', 'w') as out:
        out_t.to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    collectCoverStatistics(BAD=2)
