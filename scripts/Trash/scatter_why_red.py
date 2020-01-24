import sys
import numpy as np
import os
from scipy import stats
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


def get_color(row):
    if abs(row['log_fc']) < np.log10(fc_tr):  # or abs(row['log_pv']) < np.log10(2):
        return 'grey'
    # if row['fdrp_bh_ref'] < fdr_tr and row['fdrp_bh_alt'] < fdr_tr:
    #     return 'purple'
    if row['log_fc'] * row['log_pv'] > 0:
        return 'blue'
    else:
        return 'red'


if __name__ == '__main__':
    pt_binom = pd.read_table(os.path.expanduser("~/scatter_why_red/CTCF_new_nb.tsv"))
    pt_neg = pd.read_table(os.path.expanduser("~/scatter_why_red/CTCF_NEGBINOM.tsv"))

    perf_tr = 0.0005
    fc_tr = 4
    pvalue_tr = 1.3
    fdr_tr = 0.05

    pt_binom = pt_binom[
        (pt_binom['motif_log_pref'] >= -np.log10(perf_tr)) & (pt_binom['motif_log_palt'] >= -np.log10(perf_tr))]
    # pt_binom = pt_binom[~(pt_binom['fdrp_bh_alt'].isnull() | pt_binom['fdrp_bh_ref'].isnull())]
    pt_binom = pt_binom[(pt_binom['fdrp_bh_alt'] <= fdr_tr) | (pt_binom['fdrp_bh_ref'] <= fdr_tr)]
    pt_binom['log_pv'] = (np.log10(
        pt_binom[['fdrp_bh_ref', 'fdrp_bh_alt']]).min(axis=1)) \
                         * np.sign(pt_binom['fdrp_bh_alt'] - pt_binom['fdrp_bh_ref'])

    pt_binom['log_fc'] = pt_binom['fold_change']
    pt_binom['col'] = pt_binom.apply(lambda x: get_color(x), axis=1)

    pt_neg = pt_neg[(pt_neg['motif_log_pref'] >= -np.log10(perf_tr)) & (pt_neg['motif_log_palt'] >= -np.log10(perf_tr))]
    # pt_neg = pt_neg[~(pt_neg['fdrp_bh_alt'].isnull() | pt_neg['fdrp_bh_ref'].isnull())]
    pt_neg = pt_neg[(pt_neg['fdrp_bh_alt'] <= fdr_tr) | (pt_neg['fdrp_bh_ref'] <= fdr_tr)]
    pt_neg['log_pv'] = (np.log10(
        pt_neg[['fdrp_bh_ref', 'fdrp_bh_alt']]).min(axis=1)) \
                       * np.sign(pt_neg['fdrp_bh_alt'] - pt_neg['fdrp_bh_ref'])

    pt_neg['log_fc'] = pt_neg['fold_change']
    pt_neg['col'] = pt_neg.apply(lambda x: get_color(x), axis=1)

    # pt[pt['col'] == 'red'].to_csv(os.path.expanduser('~/TF_FC/{}_false_positives.tsv'.format(
    #     name.replace('_fc.tsv', ''))), sep='\t')

    # print(pt[pt['col'] == 'red'].groupby('BAD_mostsig_ref', as_index=False).count()[['BAD_mostsig_ref', 'col']])
    # print(pt.groupby('BAD_mostsig_ref', as_index=False).count()[['BAD_mostsig_ref', 'col']])

    pt_binom_grey = pt_binom[pt_binom['col'] == 'grey']
    pt_binom_red = pt_binom[pt_binom['col'] == 'red']
    pt_binom_blue = pt_binom[pt_binom['col'] == 'blue']

    pt_neg_grey = pt_neg[pt_neg['col'] == 'grey']
    pt_neg_red = pt_neg[pt_neg['col'] == 'red']
    pt_neg_blue = pt_neg[pt_neg['col'] == 'blue']

    print(len(pt_binom_red.index) + len(pt_binom_blue.index), len(pt_neg_red.index) + len(pt_neg_blue.index))

    neg_gr = pt_neg.groupby('ID').count()['pos']
    binom_gr = pt_binom.groupby('ID').count()['pos']

    IDs = (set(pt_binom['ID']) & set(pt_neg['ID'])) \
          - set((neg_gr[neg_gr > 1]).index) \
          - set((binom_gr[binom_gr > 1]).index)

    print(len(IDs), set((pt_binom.groupby('ID').count()['pos'] > 1).index))

    pt_binom = pt_binom[pt_binom['ID'].isin(IDs)]
    pt_neg = pt_neg[pt_neg['ID'].isin(IDs)]

    pt_binom.sort_values(by=['ID'], inplace=True)
    pt_neg.sort_values(by=['ID'], inplace=True)

    print(len(pt_binom.index), len(pt_neg.index))

    pt_binom['redness'] = -pt_binom['log_pv'] * np.sign(pt_binom['log_fc'])
    pt_neg['redness'] = -pt_neg['log_pv'] * np.sign(pt_neg['log_fc'])

    fig, ax = plt.subplots(figsize=(10, 8))
    plt.scatter(pt_binom['redness'], pt_neg['redness'])
    plt.plot([-200, 200], [-200, 200], color='black')
    plt.xlabel('new_nb_redness')
    plt.ylabel('neg_binom_redness')
    plt.grid(True)
    plt.savefig(os.path.expanduser('~/scatter_why_red/newnb-nb.png'))
