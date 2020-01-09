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
            tmp_df['alt_p'] = sum_df['p_value_alt']
            tmp_df['ref_p'] = sum_df['p_value_ref']
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


if __name__ == '__main__':
    agr_dir = os.path.expanduser('~/DATA/Processed/TF_P-values/')
    CollectRS()
    CollectPValue()
