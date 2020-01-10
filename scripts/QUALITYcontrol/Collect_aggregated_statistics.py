import sys
import os
import numpy as np
import pandas as pd
import json

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import parameters_path
from scripts.HELPERS.helpers import states


def CollectPValue(mode='by'):
    if mode == 'by':
        refname = 'fdrp_by_ref'
        altname = 'fdrp_by_alt'
    elif mode == 'bh_fisher':
        refname = 'fdrp_ref_fisher'
        altname = 'fdrp_alt_fisher'
    elif mode == 'bh':
        refname = 'fdrp_ref'
        altname = 'fdrp_alt'
    else:
        raise ValueError(mode)

    out_t = None
    for file_name in os.listdir(agr_dir):
        print(file_name)
        df = pd.read_table(agr_dir + file_name)
        if df.empty:
            continue
        df = df[df['ID'] != '.']
        sum_df = df[[refname, altname]]

        if out_t is None:
            out_t = pd.DataFrame()
            out_t['alt_p'] = sum_df[altname]
            out_t['ref_p'] = sum_df[refname]
            out_t.fillna(1, inplace=True)
            out_t = out_t.groupby(['alt_p', 'ref_p']).size().reset_index(name='counts')
        else:
            tmp_df = pd.DataFrame()
            tmp_df['alt_p'] = sum_df[altname]
            tmp_df['ref_p'] = sum_df[refname]
            tmp_df.fillna(1, inplace=True)
            tmp_df = tmp_df.groupby(['alt_p', 'ref_p']).size().reset_index(name='counts')
            out_t = out_t.append(tmp_df).groupby(['alt_p', 'ref_p'], as_index=False).sum()
            print(out_t)
    if out_t is None:
        return
    with open(parameters_path + 'fdr_pvalue_bias_statistics_{}.tsv'.format(mode), 'w') as out:
        out_t.to_csv(out, sep="\t", index=False)


def CollectRS():
    out_set = set()
    for file_name in os.listdir(agr_dir):
        print(file_name)
        df = pd.read_table(agr_dir + file_name)
        if df.empty:
            continue
        df = df[df['ID'] != '.']
        out_set |= set(df['ID'])
    with open(parameters_path + 'RS_list.tsv', 'w') as fd:
        json.dump(list(out_set), fd)


def CollectMaxCover():
    out_t = None
    for file_name in os.listdir(agr_dir):
        print(file_name)
        df = pd.read_table(agr_dir + file_name)
        if df.empty:
            continue
        df = df[df['ID'] != '.']
        df = df[(df['fdrp_by_ref'] <= 0.05) | (df['fdrp_by_alt'] <= 0.05)]
        sum_df = df[['max_cover']]

        if out_t is None:
            out_t = pd.DataFrame()
            out_t['max_cover'] = sum_df['max_cover']
            out_t.fillna(1, inplace=True)
            out_t = out_t.groupby(['max_cover']).size().reset_index(name='counts')
        else:
            tmp_df = pd.DataFrame()
            tmp_df['max_cover'] = sum_df['max_cover']
            tmp_df.fillna(1, inplace=True)
            tmp_df = tmp_df.groupby(['max_cover']).size().reset_index(name='counts')
            out_t = out_t.append(tmp_df).groupby(['max_cover'], as_index=False).sum()
            print(out_t)
    if out_t is None:
        return
    with open(parameters_path + 'fdr_mc_bias_statistics.tsv', 'w') as out:
        out_t.to_csv(out, sep="\t", index=False)


def CollectEffectSize(mode='maxdepth'):
    out_t1 = None
    out_t2 = None
    for file_name in os.listdir(agr_dir):
        print(file_name)
        df = pd.read_table(agr_dir + file_name)
        if df.empty:
            continue
        df = df[df['ID'] != '.']
        df1 = df[(df['fdrp_by_ref'] <= 0.05) | (df['fdrp_by_alt'] <= 0.05)]
        df2 = df[(df['fdrp_by_ref'] > 0.05) & (df['fdrp_by_alt'] > 0.05)]
        sum_df1 = df1[['m_' + mode]]
        sum_df2 = df2[['m_' + mode]]

        if out_t1 is None:
            out_t1 = pd.DataFrame()
            out_t1['metric'] = sum_df1['m_' + mode]
            out_t1.fillna(0, inplace=True)
            out_t1 = out_t1.groupby(['metric']).size().reset_index(name='counts')
        else:
            tmp_df = pd.DataFrame()
            tmp_df['metric'] = sum_df1['m_' + mode]
            tmp_df.fillna(0, inplace=True)
            tmp_df = tmp_df.groupby(['metric']).size().reset_index(name='counts')
            out_t1 = out_t1.append(tmp_df).groupby(['metric'], as_index=False).sum()
            print(out_t1)

        if out_t2 is None:
            out_t2 = pd.DataFrame()
            out_t2['metric'] = sum_df2['m_' + mode]
            out_t2.fillna(0, inplace=True)
            out_t2 = out_t2.groupby(['metric']).size().reset_index(name='counts')
        else:
            tmp_df = pd.DataFrame()
            tmp_df['metric'] = sum_df2['m_' + mode]
            tmp_df.fillna(0, inplace=True)
            tmp_df = tmp_df.groupby(['metric']).size().reset_index(name='counts')
            out_t2 = out_t2.append(tmp_df).groupby(['metric'], as_index=False).sum()
            print(out_t2)

    if out_t1 is None or out_t2 is None:
        return
    with open(parameters_path + 'fdr_effect_size_le_005.tsv', 'w') as out:
        out_t1.to_csv(out, sep="\t", index=False)

    with open(parameters_path + 'fdr_effect_size_gr_005.tsv', 'w') as out:
        out_t2.to_csv(out, sep="\t", index=False)


if __name__ == '__main__':
    agr_dir = os.path.expanduser('~/RESULTS/release-100120_Dipper/TF_P-values/')
    # CollectRS()
    # CollectPValue()
    # CollectMaxCover()
    CollectEffectSize()
