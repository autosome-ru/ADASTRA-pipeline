import sys
import numpy as np
import os
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


def get_color(row):
    if abs(row['log_fc']) < np.log10(fc_tr): # or abs(row['log_pv']) < np.log10(2):
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
        pt = pt[(pt['fdrp_by_alt'] <= fdr_tr) | (pt['fdrp_by_ref'] <= fdr_tr)]
        pt['log_pv'] = (np.log10(
            pt['fdrp_by_ref']) - np.log10(pt['fdrp_by_alt']))

        # pt['log_pv'] = (- np.log10(pt['fdrp_by_alt']))

        # pt['log_pv'] = (np.log10(
        #     pt[['fdrp_by_ref', 'fdrp_by_alt']]).min(axis=1)) \
        #                * np.sign(pt['fdrp_by_alt'] - pt['fdrp_by_ref'])
        pt['log_fc'] = pt['fold_change']
        pt['col'] = pt.apply(lambda x: get_color(x), axis=1)

        blue = len(pt[pt['col'] == "blue"].index)
        red = len(pt[pt['col'] == "red"].index)
        grey = len(pt[pt['col'] == "grey"].index)

        pt[pt['col'] == 'red'].to_csv(os.path.expanduser('~/TF_FC/{}_false_positives.tsv'.format(
            name.replace('_fc.tsv', ''))), sep='\t')

        print(pt[pt['col'] == 'red'].groupby('BAD_mostsig_ref', as_index=False).count()[['BAD_mostsig_ref', 'col']])
        print(pt.groupby('BAD_mostsig_ref', as_index=False).count()[['BAD_mostsig_ref', 'col']])

        fig, ax = plt.subplots(figsize=(10, 10))
        pt_grey = pt[pt['col'] == 'grey']
        pt_red = pt[pt['col'] == 'red']
        pt_blue = pt[pt['col'] == 'blue']
        sns.scatterplot(x=pt_grey['log_pv'], y=pt_grey['log_fc'], color='grey')
        sns.scatterplot(x=pt_red['log_pv'], y=pt_red['log_fc'], color='red')
        sns.scatterplot(x=pt_blue['log_pv'], y=pt_blue['log_fc'], color='blue')
        plt.grid(True)
        plt.title(name + '\n' +
                  '\n'.join(['perfectos_pv <= {}', 'fold_change >= {}', 'fdr <= {}', fix]).format(perf_tr, fc_tr,
                                                                                                  fdr_tr))
        plt.xlabel('delta -log10 fdr_p')
        plt.ylabel('log10 prefectos foldchange')
        label = 'blue/red: {}/{}({:.1f}%),\ngrey/all={:.1f}%'.format(blue, red,
                                                                     100 * blue / (blue + red),
                                                                     100 * grey / (grey + blue + red))
        plt.text(x=max(pt['log_pv']) / 5, y=min(pt['log_fc']) / 2, s=label)
        plt.text(x=min(pt['log_pv']) * 4 / 5, y=min(pt['log_fc']) / 2,
                 s=' P-value treshold: {},\nlog FC treshold: {}'.format(round(fdr_tr, 2), round(fc_tr, 1)))

        plt.savefig(os.path.expanduser("~/TF_FC/{}_p_tr={:.2f}_fc_tr={:.2f}_fdr_tr={:.2f}.png".format(
            name.replace('_fc.tsv', ''), perf_tr, fc_tr, fdr_tr)))
