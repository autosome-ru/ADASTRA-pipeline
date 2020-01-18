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

        blue_fr_ref = []
        blue_ct_ref = []
        red_ct_ref = []

        blue_fr_alt = []
        blue_ct_alt = []
        red_ct_alt = []

        x = np.linspace(0.001, 0.1, 100)

        for fdr_tr in x:
            print(fdr_tr)

            # pt['log_pv'] = (- np.log10(pt['fdrp_by_alt']))

            pt2 = pt[(pt['fdrp_by_alt'] <= fdr_tr)]

            blue_alt = len(pt2[(pt2['col'] == "blue")].index)
            red_alt = len(pt2[pt2['col'] == "red"].index)
            grey_alt = len(pt2[pt2['col'] == "grey"].index)

            pt2 = pt[(pt['fdrp_by_ref'] <= fdr_tr)]

            blue_ref = len(pt2[(pt2['col'] == "blue")].index)
            red_ref = len(pt2[pt2['col'] == "red"].index)
            grey_ref = len(pt2[pt2['col'] == "grey"].index)

            blue_fr_ref.append(blue_ref / (blue_ref + red_ref))
            blue_ct_ref.append(blue_ref)
            red_ct_ref.append(red_ref)

            blue_fr_alt.append(blue_alt / (blue_alt + red_alt))
            blue_ct_alt.append(blue_alt)
            red_ct_alt.append(red_alt)

        fig, ax = plt.subplots(figsize=(10, 8))
        plt.scatter(x, blue_fr_ref, color='C0', label='ref')
        plt.scatter(x, blue_fr_alt, color='C1', label='alt')
        plt.grid(True)
        plt.legend()
        plt.xlabel('fdr_tr')
        plt.ylabel('blue / (blue + red)')
        plt.title('CTCF blue fraction')
        plt.savefig(os.path.expanduser("~/TF_FC/{}blue_frac_refalt_p_tr={:.2f}_fc_tr={:.2f}.png".format(
            name.replace('_fc.tsv', ''), perf_tr, fc_tr, fdr_tr)))

        fig, ax = plt.subplots(figsize=(10, 8))
        plt.scatter(x, blue_ct_ref, color='C0', label='ref')
        plt.scatter(x, blue_ct_alt, color='C1', label='alt')
        plt.grid(True)
        plt.legend()
        plt.xlabel('fdr_tr')
        plt.ylabel('count')
        plt.title('CTCF blue counts')
        plt.savefig(os.path.expanduser("~/TF_FC/{}blue_counts_refalt_p_tr={:.2f}_fc_tr={:.2f}.png".format(
            name.replace('_fc.tsv', ''), perf_tr, fc_tr, fdr_tr)))

        fig, ax = plt.subplots(figsize=(10, 8))
        plt.scatter(x, red_ct_ref, color='C0', label='ref')
        plt.scatter(x, red_ct_alt, color='C1', label='alt')
        plt.grid(True)
        plt.legend()
        plt.xlabel('fdr_tr')
        plt.ylabel('count')
        plt.title('CTCF red counts')
        plt.savefig(os.path.expanduser("~/TF_FC/{}red_counts_refalt_p_tr={:.2f}_fc_tr={:.2f}.png".format(
            name.replace('_fc.tsv', ''), perf_tr, fc_tr, fdr_tr)))



            # fig, ax = plt.subplots(figsize=(10, 10))
            # pt_grey = pt[pt['col'] == 'grey']
            # pt_red = pt[pt['col'] == 'red']
            # pt_blue = pt[pt['col'] == 'blue']
            # sns.scatterplot(x=pt_grey['log_pv'], y=pt_grey['log_fc'], color='grey')
            # sns.scatterplot(x=pt_red['log_pv'], y=pt_red['log_fc'], color='red')
            # sns.scatterplot(x=pt_blue['log_pv'], y=pt_blue['log_fc'], color='blue')
            # plt.grid(True)
            # plt.title(name + '\n' +
            #           '\n'.join(['perfectos_pv <= {}', 'fold_change >= {}', 'fdr <= {}', fix]).format(perf_tr, fc_tr,
            #                                                                                           fdr_tr))
            # plt.xlabel('delta -log10 fdr_p')
            # plt.ylabel('log10 prefectos foldchange')
            # label = 'blue/red: {}/{}({:.1f}%),\ngrey/all={:.1f}%'.format(blue, red,
            #                                                              100 * blue / (blue + red),
            #                                                              100 * grey / (grey + blue + red))
            # plt.text(x=max(pt['log_pv']) / 5, y=min(pt['log_fc']) / 2, s=label)
            # plt.text(x=min(pt['log_pv']) * 4 / 5, y=min(pt['log_fc']) / 2,
            #          s=' P-value treshold: {},\nlog FC treshold: {}'.format(round(fdr_tr, 2), round(fc_tr, 1)))
            #
            # plt.savefig(os.path.expanduser("~/TF_FC/{}_p_tr={:.2f}_fc_tr={:.2f}_fdr_tr={:.2f}.png".format(
            #     name.replace('_fc.tsv', ''), perf_tr, fc_tr, fdr_tr)))
