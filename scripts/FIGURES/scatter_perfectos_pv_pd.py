import sys
import numpy as np
import os
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


def get_color(row):
    if abs(row['log_fc']) < np.log10(fc_tr):  # or abs(row['log_pv']) < np.log10(2):
        return 'grey'
    if row['log_fc'] * row['log_pv'] > 0:
        return 'blue'
    else:
        return 'red'


if __name__ == '__main__':
    top10_names = ['CTCF_HUMAN_fc.tsv', 'FOXA1_HUMAN_fc.tsv', 'ESR1_HUMAN_fc.tsv', 'SPI1_HUMAN_fc.tsv',
                   'ANDR_HUMAN_fc.tsv', 'STAT1_HUMAN_fc.tsv', 'IKZF1_HUMAN_fc.tsv', 'CEBPB_HUMAN_fc.tsv',
                   'ZFX_HUMAN_fc.tsv', 'CREB1_HUMAN_fc.tsv']
    for name in ['CTCF_HUMAN_fc.tsv']:
        pt = pd.read_table(os.path.expanduser("~/TF_FC/{}".format(name)))

        perf_tr = 0.0005
        fc_tr = 4
        pvalue_tr = 1.3
        fdr_tr = 0.05
        fix = ''

        pt = pt[(pt['motif_log_pref'] >= -np.log10(perf_tr)) & (pt['motif_log_palt'] >= -np.log10(perf_tr))]
        pt = pt[~(pt['fdrp_by_alt'].isnull() | pt['fdrp_by_ref'].isnull())]
        # pt['log_pv'] = (np.log10(
        #     pt['fdrp_by_ref']) - np.log10(pt['fdrp_by_alt']))

        pt['log_pv'] = (np.log10(
            pt[['fdrp_by_ref', 'fdrp_by_alt']]).min(axis=1)) \
                       * np.sign(pt['fdrp_by_alt'] - pt['fdrp_by_ref'])
        pt['log_fc'] = pt['fold_change']
        pt['col'] = pt.apply(lambda x: get_color(x), axis=1)

        blue_fr = []
        blue_ct = []
        red_ct = []

        #x = np.linspace(0.001, 0.1, 100)
        x = np.linspace(30, 1500, 250)

        #for fdr_tr in x:
        for maxcov in x:
            # pt['log_pv'] = (- np.log10(pt['fdrp_by_alt']))

            pt = pt[(pt['fdrp_by_alt'] <= fdr_tr) | (pt['fdrp_by_ref'] <= fdr_tr)]
            #pt2 = pt[(pt['altc_mostsig_alt'] <= maxcov) | (pt['refc_mostsig_ref'] <= maxcov)]
            pt2 = pt[(pt['max_cover'] <= maxcov)]

            blue = len(pt2[(pt2['col'] == "blue")].index)
            red = len(pt2[pt2['col'] == "red"].index)
            grey = len(pt2[pt2['col'] == "grey"].index)

            blue_fr.append(blue / (blue + red))
            blue_ct.append(blue)
            red_ct.append(red)

        fig, ax = plt.subplots(figsize=(10, 8))
        plt.scatter(x, blue_fr, color='black')
        plt.grid(True)
        # plt.xlabel('fdr_tr')
        plt.xlabel('max cov tr')
        plt.ylabel('blue / (blue + red)')
        plt.title('CTCF blue fraction')
        plt.savefig(os.path.expanduser("~/TF_FC/{}blue_frac_p_tr={:.2f}_fc_tr={:.2f}.png".format(
            name.replace('_fc.tsv', ''), perf_tr, fc_tr, fdr_tr)))

        fig, ax = plt.subplots(figsize=(10, 8))
        plt.scatter(x, blue_ct, color='blue')
        plt.scatter(x, red_ct, color='red')
        plt.grid(True)
        # plt.xlabel('fdr_tr')
        plt.xlabel('max cov tr')
        plt.ylabel('count')
        plt.title('CTCF blue red counts')
        plt.savefig(os.path.expanduser("~/TF_FC/{}blue_red_counts_refalt_p_tr={:.2f}_fc_tr={:.2f}.png".format(
            name.replace('_fc.tsv', ''), perf_tr, fc_tr, fdr_tr)))
