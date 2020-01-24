import sys
import numpy as np
import os
from scipy import stats
import pandas as pd
from matplotlib import pyplot as plt
sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.FIGURES import style_config


def get_color(row):
    if abs(row['log_fc']) < np.log10(fc_tr) or abs(row['log_pv']) < -np.log10(fdr_tr):
        return '#CCCCCC'
    # if row[field + '_ref'] < fdr_tr and row[field + '_alt'] < fdr_tr:
    #     return 'purple'
    if row['log_fc'] * row['log_pv'] > 0:
        return '#005AB5'
    else:
        return '#DC3220'


if __name__ == '__main__':
    top10_names = ['CTCF_HUMAN.tsv', 'SPI1_HUMAN.tsv',
                   'FOXA1_HUMAN.tsv', 'CEBPB_HUMAN.tsv',
                   'ANDR_HUMAN.tsv', 'ESR1_HUMAN.tsv',
                   'NRF1_HUMAN.tsv', 'DUX4_HUMAN.tsv',
                   'CREB1_HUMAN.tsv', 'AP2A_HUMAN.tsv']
    for name in top10_names:
        pt = pd.read_table(os.path.expanduser("~/scatter_why_red/{}".format(name)))

        field = 'fdrp_bh'

        perf_tr = 0.0005
        fc_tr = 2
        fdr_tr = 0.05

        ##
        # pt = pt[pt['BAD'] != 4/3]


        pt = pt[(pt['motif_log_pref'] >= -np.log10(perf_tr)) & (pt['motif_log_palt'] >= -np.log10(perf_tr))]
        pt = pt[~(pt[field + '_alt'].isnull() | pt[field + '_ref'].isnull())]
        #pt = pt[(pt[field + '_alt'] <= fdr_tr) | (pt[field + '_ref'] <= fdr_tr)]


        # pt['log_pv'] = (np.log10(
        #     pt[field + '_ref']) - np.log10(pt[field + '_alt']))

        # pt['log_pv'] = (- np.log10(pt[field + '_alt']))

        # pt['log_pv'] = (np.log10(
        #     pt[[field + '_ref', field + '_alt']]).min(axis=1)) \
        #                * np.sign(pt[field + '_alt'] - pt[field + '_ref'])
        #
        pt['log_pv'] = (np.log10(
            pt[[field + '_ref', field + '_alt']]).min(axis=1)) \
                       * np.sign(pt[field + '_alt'] - pt[field + '_ref'])


        pt['log_fc'] = pt['fold_change']
        pt['col'] = pt.apply(lambda x: get_color(x), axis=1)
        #pt = pt[pt['max_cover'] <= maxcov_tr]



        #additional columns
        # pt['minc_mostsig_ref'] = pt[['refc_mostsig_ref', 'altc_mostsig_ref']].min(axis=1)
        # pt['minc_mostsig_alt'] = pt[['refc_mostsig_alt', 'altc_mostsig_alt']].min(axis=1)
        # pt['maxc_mostsig_ref'] = pt[['refc_mostsig_ref', 'altc_mostsig_ref']].max(axis=1)
        # pt['maxc_mostsig_alt'] = pt[['refc_mostsig_alt', 'altc_mostsig_alt']].max(axis=1)
        # pt[field + '_min'] = pt[[field + '_ref', field + '_alt']].min(axis=1)

        blue = len(pt[pt['col'] == "#005AB5"].index)
        red = len(pt[pt['col'] == "#DC3220"].index)
        grey = len(pt[pt['col'] == "#CCCCCC"].index)

        # i = len(pt[(pt['col'] == "blue") & (pt['log_pv'] > 0)].index)
        # ii = len(pt[(pt['col'] == "red") & (pt['log_pv'] < 0)].index)
        # iii = len(pt[(pt['col'] == "blue") & (pt['log_pv'] < 0)].index)
        # iv = len(pt[(pt['col'] == "red") & (pt['log_pv'] > 0)].index)
        #
        # print(i, ii, iii, iv)

        pt[pt['col'] == "#DC3220"].to_csv(os.path.expanduser('~/TF_FC/{}_false_positives.tsv'.format(
            name.replace('_fc.tsv', ''))), sep='\t')

        # print(pt[pt['col'] == 'red'].groupby('BAD_mostsig_ref', as_index=False).count()[['BAD_mostsig_ref', 'col']])
        # print(pt.groupby('BAD_mostsig_ref', as_index=False).count()[['BAD_mostsig_ref', 'col']])


        pt_grey = pt[pt['col'] == "#CCCCCC"]
        pt_red = pt[pt['col'] == "#DC3220"]
        pt_blue = pt[pt['col'] == "#005AB5"]
        # pt_purple = pt[pt['col'] == 'purple']

        # fig, ax = plt.subplots(figsize=(10, 8))
        #
        # # print(pt_red.groupby('max_cover', as_index=False).count()[['max_cover', 'col']])
        # # print(pt_blue.groupby('max_cover', as_index=False).count()[['max_cover', 'col']])
        #
        #
        #
        # plt.hist(pt['max_cover'], color='black', alpha=0.5, bins=45, range=(10, 100))
        # # plt.hist(pt_red['max_cover'], color='red', alpha=0.5)
        # plt.grid(True)
        # plt.xlabel('max_cover')
        # plt.ylabel('count')
        # plt.title(name + ' ' + field + '\n' +
        #           '\n'.join(['perfectos_pv <= {}', 'fold_change >= {}', 'fdr <= {}']).format(perf_tr, fc_tr,
        #                                                                                      fdr_tr))
        # plt.show()

        # print('\n'.join(map(str, list(pt_red['ID']))))
        # print('blue')
        # print('\n'.join(map(str, list(pt_blue['ID']))))

        # print(pt.info())
        # wc = {}
        # for col in pt.columns:
        #     if pt[col].dtype == np.float_ or pt[col].dtype == np.int_:
        #         wc[col] = (stats.ranksums(pt_blue[col], pt_red[col])[1],
        #                    np.nanmedian(pt_blue[col]), np.nanmedian(pt_red[col]))
        # l = list(wc.items())
        # l = sorted(l, key=lambda x: x[1])
        # print('\n'.join(map(str, l)))
        #
        # print(pt[pt['ID'] == '.'])
        fig, ax = plt.subplots(figsize=(10, 8))
        plt.scatter(x=pt['log_pv'], y=pt['log_fc'], c=pt['col'], s=5)
        # sns.scatterplot(x=pt_grey['log_pv'], y=pt_grey['log_fc'], color='grey')
        # sns.scatterplot(x=pt_red['log_pv'], y=pt_red['log_fc'], color='red')
        # sns.scatterplot(x=pt_blue['log_pv'], y=pt_blue['log_fc'], color='blue')
        plt.grid(True)
        plt.title(name + ' ' +field + '\n' +
                  '\n'.join(['perfectos_pv <= {}', 'fold_change >= {}', 'fdr <= {}']).format(perf_tr, fc_tr,
                                                                                                  fdr_tr))
        plt.xlabel('best -log10 fdr_p')
        plt.ylabel('log10 prefectos foldchange')
        label = 'blue/red: {}/{}({:.1f}%),\ngrey/all={:.1f}%'.format(blue, red,
                                                                     100 * blue / (blue + red),
                                                                     100 * grey / (grey + blue + red))
        plt.text(x=max(pt['log_pv']) / 5, y=min(pt['log_fc']) / 2, s=label)
        plt.text(x=min(pt['log_pv']) * 4 / 5, y=min(pt['log_fc']) / 2,
                 s=' P-value treshold: {},\nFC treshold: {}'.format(round(fdr_tr, 2), round(fc_tr, 1)))

        plt.savefig(os.path.expanduser("~/TF_FC/final/{}_p_tr={:.2f}_fc_tr={:.2f}_fdr_tr={:.2f}_{}.png".format(
            name.replace('_fc.tsv', ''), perf_tr, fc_tr, fdr_tr, field)))

