import sys
import numpy as np
import os
from scipy import stats
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
        maxcov_tr = 500000
        fix = ''

        pt = pt[(pt['motif_log_pref'] >= -np.log10(perf_tr)) & (pt['motif_log_palt'] >= -np.log10(perf_tr))]
        pt = pt[~(pt['fdrp_by_alt'].isnull() | pt['fdrp_by_ref'].isnull())]
        pt = pt[(pt['fdrp_by_alt'] <= fdr_tr) | (pt['fdrp_by_ref'] <= fdr_tr)]


        # pt['log_pv'] = (np.log10(
        #     pt['fdrp_by_ref']) - np.log10(pt['fdrp_by_alt']))

        # pt['log_pv'] = (- np.log10(pt['fdrp_by_alt']))

        pt['log_pv'] = (np.log10(
            pt[['fdrp_by_ref', 'fdrp_by_alt']]).min(axis=1)) \
                       * np.sign(pt['fdrp_by_alt'] - pt['fdrp_by_ref'])
        pt['log_fc'] = pt['fold_change']
        pt['col'] = pt.apply(lambda x: get_color(x), axis=1)
        pt = pt[pt['max_cover'] <= maxcov_tr]


        #additional columns
        pt['minc_mostsig_ref'] = pt[['refc_mostsig_ref', 'altc_mostsig_ref']].min(axis=1)
        pt['minc_mostsig_alt'] = pt[['refc_mostsig_alt', 'altc_mostsig_alt']].min(axis=1)
        pt['maxc_mostsig_ref'] = pt[['refc_mostsig_ref', 'altc_mostsig_ref']].max(axis=1)
        pt['maxc_mostsig_alt'] = pt[['refc_mostsig_alt', 'altc_mostsig_alt']].max(axis=1)
        pt['logitp_min'] = pt[['logitp_ref', 'logitp_alt']].min(axis=1)

        blue = len(pt[pt['col'] == "blue"].index)
        red = len(pt[pt['col'] == "red"].index)
        grey = len(pt[pt['col'] == "grey"].index)

        i = len(pt[(pt['col'] == "blue") & (pt['log_pv'] > 0)].index)
        ii = len(pt[(pt['col'] == "red") & (pt['log_pv'] < 0)].index)
        iii = len(pt[(pt['col'] == "blue") & (pt['log_pv'] < 0)].index)
        iv = len(pt[(pt['col'] == "red") & (pt['log_pv'] > 0)].index)

        print(i, ii, iii, iv)

        pt[pt['col'] == 'red'].to_csv(os.path.expanduser('~/TF_FC/{}_false_positives.tsv'.format(
            name.replace('_fc.tsv', ''))), sep='\t')

        # print(pt[pt['col'] == 'red'].groupby('BAD_mostsig_ref', as_index=False).count()[['BAD_mostsig_ref', 'col']])
        # print(pt.groupby('BAD_mostsig_ref', as_index=False).count()[['BAD_mostsig_ref', 'col']])

        fig, ax = plt.subplots(figsize=(10, 8))
        pt_grey = pt[pt['col'] == 'grey']
        pt_red = pt[pt['col'] == 'red']
        pt_blue = pt[pt['col'] == 'blue']

        print('\n'.join(map(str, list(pt_red['ID']))))
        print('blue')
        print('\n'.join(map(str, list(pt_blue['ID']))))


        # print(pt.info())
        wc = {}
        for col in pt.columns:
            if pt[col].dtype == np.float_ or pt[col].dtype == np.int_:
                wc[col] = (stats.ranksums(pt_blue[col], pt_red[col])[1],
                           np.nanmedian(pt_blue[col]), np.nanmedian(pt_red[col]))
        l = list(wc.items())
        l = sorted(l, key=lambda x: x[1])
        print('\n'.join(map(str, l)))

        plt.scatter(x=pt['log_pv'], y=pt['log_fc'], c=pt['col'], s=5)
        # sns.scatterplot(x=pt_grey['log_pv'], y=pt_grey['log_fc'], color='grey')
        # sns.scatterplot(x=pt_red['log_pv'], y=pt_red['log_fc'], color='red')
        # sns.scatterplot(x=pt_blue['log_pv'], y=pt_blue['log_fc'], color='blue')
        plt.grid(True)
        plt.title(name + '\n' +
                  '\n'.join(['perfectos_pv <= {}', 'fold_change >= {}', 'fdr <= {}', 'max_cov <= {}']).format(perf_tr, fc_tr,
                                                                                                  fdr_tr, maxcov_tr))
        plt.xlabel('delta -log10 fdr_p')
        plt.ylabel('log10 prefectos foldchange')
        label = 'blue/red: {}/{}({:.1f}%),\ngrey/all={:.1f}%'.format(blue, red,
                                                                     100 * blue / (blue + red),
                                                                     100 * grey / (grey + blue + red))
        plt.text(x=max(pt['log_pv']) / 5, y=min(pt['log_fc']) / 2, s=label)
        plt.text(x=min(pt['log_pv']) * 4 / 5, y=min(pt['log_fc']) / 2,
                 s=' P-value treshold: {},\nFC treshold: {}'.format(round(fdr_tr, 2), round(fc_tr, 1)))

        plt.savefig(os.path.expanduser("~/TF_FC/{}_p_tr={:.2f}_fc_tr={:.2f}_fdr_tr={:.2f}_maxcov_tr={}.png".format(
            name.replace('_fc.tsv', ''), perf_tr, fc_tr, fdr_tr, maxcov_tr)))
