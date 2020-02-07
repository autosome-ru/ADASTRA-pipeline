import sys
import os.path
import json
import pandas as pd
import numpy as np

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths_for_components import parameters_path, cl_dict_path


def collectRefAltStatistics(key_name=None, BAD=None):
    with open(cl_dict_path, "r") as read_file:
        cell_lines_dict = json.loads(read_file.readline())
    out_t = None

    for key in cell_lines_dict:
        if key_name is not None:
            if key not in key_name:  # <------
                continue
        for align_path in cell_lines_dict[key]:
            align_path = align_path.replace('table_p', 'table_BADs')
            if not os.path.isfile(align_path):
                continue
            print(align_path)
            df = pd.read_table(align_path)
            if df.empty:
                continue
            df = df[df['ID'] != '.']
            if BAD is not None:
                sum_df = df[df['BAD'] == BAD][['ref_read_counts', 'alt_read_counts']]  # <------
            else:
                sum_df = df[['ref_read_counts', 'alt_read_counts']]

            if out_t is None:
                out_t = pd.DataFrame()
                s1 = sum_df['ref_read_counts'].value_counts()
                s2 = sum_df['alt_read_counts'].value_counts()
                out_t['allele_reads'] = sorted(list(set(s1.index) | set(s2.index)))
                out_t.index = out_t['allele_reads']
                out_t['ref'] = s1
                out_t['alt'] = s2
                out_t = out_t.reset_index(drop=True)
                out_t.fillna(0, inplace=True)
                out_t = out_t.astype(int)
            else:
                tmp_df = pd.DataFrame()
                s1 = sum_df['ref_read_counts'].value_counts()
                s2 = sum_df['alt_read_counts'].value_counts()
                tmp_df['allele_reads'] = sorted(list(set(s1.index) | set(s2.index)))
                tmp_df.index = tmp_df['allele_reads']
                tmp_df['ref'] = s1
                tmp_df['alt'] = s2
                tmp_df = tmp_df.reset_index(drop=True)
                tmp_df.fillna(0, inplace=True)
                tmp_df = tmp_df.astype(np.int_)
                out_t = out_t.append(tmp_df).groupby('allele_reads', as_index=False).sum()
    if out_t is None:
        return
    with open(parameters_path + 'bias_statistics_BAD={:.1f}.tsv'.format(BAD), 'w') as out:
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
                sum_df = df[df['BAD'] == BAD][['ref_read_counts', 'alt_read_counts']]  # <------
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
    with open(parameters_path + 'cover_bias_statistics_BAD={:.1f}.tsv'.format(BAD), 'w') as out:
        out_t.to_csv(out, sep="\t", index=False)


def collectNegativeBinomStatistics(key_name=None, BAD=None, alt=True):
    with open(cl_dict_path, "r") as read_file:
        cell_lines_dict = json.loads(read_file.readline())
    out_t = None
    if alt:
        alt = "alt"
    else:
        alt = "ref"

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
            if BAD is not None:
                sum_df = df[df['BAD'] == BAD][['{}_read_counts'.format(alt)]]  # <------
            else:
                sum_df = df[['{}_read_counts'.format(alt)]]

            if out_t is None:
                out_t = pd.DataFrame()
                out_t['{}_counts'.format(alt)] = sum_df['{}_read_counts'.format(alt)]
                out_t = out_t.groupby(['{}_counts'.format(alt)]).size().reset_index(name='counts')
                out_t.fillna(0, inplace=True)
            else:
                tmp_df = pd.DataFrame()
                tmp_df['{}_counts'.format(alt)] = sum_df['{}_read_counts'.format(alt)]
                tmp_df = tmp_df.groupby(['{}_counts'.format(alt)]).size().reset_index(name='counts')
                tmp_df.fillna(0, inplace=True)
                out_t = out_t.append(tmp_df).groupby(['{}_counts'.format(alt)], as_index=False).sum()
    if out_t is None:
        return
    with open(parameters_path + '{}_bias_statistics_BAD={:.1f}.tsv'.format(alt, BAD), 'w') as out:
        out_t.to_csv(out, sep="\t", index=False)


def collectFixedAltStatistics(key_name=None, BAD=None):
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
            df = df[df['ID'] != '.']
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
    with open(parameters_path + 'fixed_alt_bias_statistics_BAD={:.1f}.tsv'.format(BAD), 'w') as out:
        out_t.to_csv(out, sep="\t", index=False)


def collectFixedMinStatistics(key_name=None, BAD=None):
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
            df = df[df['ID'] != '.']
            if BAD is not None:
                sum_df = df[df['BAD'] == BAD][['ref_read_counts', 'alt_read_counts']]  # <------
            else:
                sum_df = df[['ref_read_counts', 'alt_read_counts']]

            sum_df['min_read_counts'] = sum_df.min(axis=1)
            sum_df['max_read_counts'] = sum_df.max(axis=1)

            if out_t is None:
                out_t = pd.DataFrame()
                out_t['min_counts'] = sum_df['min_read_counts']
                out_t['max_counts'] = sum_df['max_read_counts']
                out_t = out_t.groupby(['min_counts', 'max_counts']).size().reset_index(name='counts')
                out_t.fillna(0, inplace=True)
            else:
                tmp_df = pd.DataFrame()
                tmp_df['min_counts'] = sum_df['min_read_counts']
                tmp_df['max_counts'] = sum_df['max_read_counts']
                tmp_df = tmp_df.groupby(['min_counts', 'max_counts']).size().reset_index(name='counts')
                tmp_df.fillna(0, inplace=True)
                out_t = out_t.append(tmp_df).groupby(['min_counts', 'max_counts'], as_index=False).sum()
    if out_t is None:
        return
    with open(parameters_path + 'fixed_min_bias_statistics_BAD={:.1f}.tsv'.format(BAD), 'w') as out:
        out_t.to_csv(out, sep="\t", index=False)


def collectPValueStatistics(key_name=None, BAD=None):
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
            df = df[df['ID'] != '.']

            assert len(df[(df['p_value_ref'] == 0) | (df['p_value_alt'] == 0)].index) == 0

            if BAD is not None:
                sum_df = df[df['BAD'] == BAD][['p_value_ref', 'p_value_alt']]  # <------
            else:
                sum_df = df[['p_value_ref', 'p_value_alt']]

            if out_t is None:
                out_t = pd.DataFrame()
                out_t['alt_p'] = sum_df['p_value_alt']
                out_t['ref_p'] = sum_df['p_value_ref']
                out_t.fillna(1, inplace=True)
                out_t = out_t.groupby(['alt_p', 'ref_p']).size().reset_index(name='counts')
            else:
                tmp_df = pd.DataFrame()
                tmp_df['alt_p'] = sum_df['p_value_alt']
                tmp_df['ref_p'] = sum_df['p_value_ref']
                tmp_df.fillna(1, inplace=True)
                tmp_df = tmp_df.groupby(['alt_p', 'ref_p']).size().reset_index(name='counts')
                out_t = out_t.append(tmp_df).groupby(['alt_p', 'ref_p'], as_index=False).sum()
                print(out_t)
    if out_t is None:
        return
    with open(parameters_path + 'pvalue_bias_statistics_BAD={:.1f}.tsv'.format(BAD), 'w') as out:
        out_t.to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    for alt in {True, False}:
        for BAD in [1, 2, 3, 4, 5, 6, 4/3, 5/2, 3/2]:
            collectNegativeBinomStatistics(alt=alt, BAD=BAD)
