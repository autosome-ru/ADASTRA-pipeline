import sys
import numpy as np
import os
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


def get_color(row):
    if abs(row['log_fc']) < np.log10(fc_tr) or abs(row['log_pv']) < -np.log10(fdr_tr):
        return grey_color
    # if row[field + '_ref'] < fdr_tr and row[field + '_alt'] < fdr_tr:
    #     return 'purple'
    if row['log_fc'] * row['log_pv'] > 0:
        return blue_color
    else:
        return red_color


if __name__ == '__main__':

    top10_names = ['CTCF_HUMAN.tsv', 'SPI1_HUMAN.tsv',
                   'FOXA1_HUMAN.tsv', 'CEBPB_HUMAN.tsv',
                   'ANDR_HUMAN.tsv', 'ESR1_HUMAN.tsv',
                   'NRF1_HUMAN.tsv', 'DUX4_HUMAN.tsv',
                   'CREB1_HUMAN.tsv', 'AP2A_HUMAN.tsv']  # , 'CTCF_Mathelier.tsv']

    sns.set(font_scale=1.4, style="ticks", font="lato", palette=('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2',
                                                                 '#D55E00', '#CC79A7'))
    sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})
    plt.rcParams['font.weight'] = "medium"
    plt.rcParams['axes.labelweight'] = 'medium'
    plt.rcParams['figure.titleweight'] = 'medium'
    plt.rcParams['axes.titleweight'] = 'medium'
    plt.rcParams['figure.figsize'] = 6, 5
    plt.rcParams["legend.framealpha"] = 1
    plt.rcParams['axes.xmargin'] = 0
    plt.rcParams['axes.ymargin'] = 0
    plt.rcParams["legend.framealpha"] = 1

    field = 'fdrp_bh'

    perf_tr = 0.0005
    fc_tr = 4
    fdr_tr = 0.005

    blue_color = '#005AB5'  # 1B7837'
    red_color = '#DC3220'  # 762A83'
    grey_color = '#CCCCCC'
    bar_width = 0.8
    bar_alpha = 0.5

    # Barplot
    df = pd.read_table(os.path.expanduser("~/blue_red_stats.tsv"))
    fig, ax = plt.subplots()
    plt.tight_layout(pad=2.5)
    df['sum'] = df["red"] + df["blue"]
    df['part'] = df["blue"] / df['sum']
    df = df[df['sum'] >= 50].sort_values("sum")
    x, blue, red, blue_n, ticks = range(len(df.index)), \
                                   df["part"].tolist(), (df["red"] / df['sum']).tolist(), df["blue"].tolist(), \
                                   df["name"].apply(lambda x: x.replace("_HUMAN", ""))

    ax.barh(x, red, left=blue, color="C5", height=bar_width, linewidth=0, tick_label=ticks, alpha=bar_alpha)
    ax.barh(x, blue, color="C4", height=bar_width, linewidth=0, alpha=bar_alpha)
    for i in x:
        ax.text(x=0.5, y=i-0.1, s='{}/{:.0f}'.format(blue_n[i], blue_n[i]/blue[i]), va="center", ha="center",
                fontdict={"size": 13})
    ax.set_xlabel('TFs sorted by # of snps in motifs')
    ax.set_ylabel('# of snps')
    ax.legend(loc='upper right')
    plt.savefig(os.path.expanduser("~/AC_7/AS_Figure_11.png"), dpi=300)
    plt.savefig(os.path.expanduser("~/AC_7/AS_Figure_11.svg"), dpi=300)

    plt.close(fig)

    # Scatters
    for name in top10_names:
        pt = pd.read_table(os.path.expanduser("~/scatter_why_red/{}".format(name)))

        pt = pt[(pt['motif_log_pref'] >= -np.log10(perf_tr)) & (pt['motif_log_palt'] >= -np.log10(perf_tr))]
        pt = pt[~(pt[field + '_alt'].isnull() | pt[field + '_ref'].isnull())]

        pt['log_pv'] = (np.log10(
            pt[[field + '_ref', field + '_alt']]).min(axis=1)) \
                       * np.sign(pt[field + '_alt'] - pt[field + '_ref'])

        pt['log_fc'] = pt['fold_change']
        pt['col'] = pt.apply(lambda x: get_color(x), axis=1)

        blue = len(pt[pt['col'] == blue_color].index)
        red = len(pt[pt['col'] == red_color].index)
        grey = len(pt[pt['col'] == grey_color].index)

        pt_grey = pt[pt['col'] == grey_color]
        pt_red = pt[pt['col'] == red_color]
        pt_blue = pt[pt['col'] == blue_color]

        fig, ax = plt.subplots()
        plt.tight_layout(pad=1.5)
        ax.scatter(x=pt['log_pv'], y=pt['log_fc'], c=pt['col'], s=5)
        ax.grid(True)
        ax.set_xlabel('Best -log10 fdr_p')
        ax.set_ylabel('Log10 motif foldchange')
        label = 'blue/red: {}/{}({:.1f}%),\ngrey/all={:.1f}%'.format(blue, red,
                                                                     100 * blue / (blue + red),
                                                                     100 * grey / (grey + blue + red))
        # plt.text(x=max(pt['log_pv']) / 5, y=min(pt['log_fc']) / 2, s=label)
        # plt.text(x=min(pt['log_pv']) * 4 / 5, y=min(pt['log_fc']) / 2,
        #          s=' P-value treshold: {},\nFC treshold: {}'.format(round(fdr_tr, 2), round(fc_tr, 1)))
        # plt.text()

        plt.savefig(os.path.expanduser("~/AC_7/AS_Figure_7_{}_p_tr={:.2f}_fc_tr={:.2f}_fdr_tr={:.2f}_{}.svg".format(
            name.replace('_fc.tsv', ''), perf_tr, fc_tr, fdr_tr, field)), dpi=300)
        plt.close(fig)
