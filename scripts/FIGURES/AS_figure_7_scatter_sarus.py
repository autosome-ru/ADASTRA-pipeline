import sys
import numpy as np
import os
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


def get_color(row):
    if abs(row['log_fc']) < fc_tr or abs(row['log_pv']) < -np.log10(fdr_tr):
        return grey_color
    # if row[field + '_ref'] < fdr_tr and row[field + '_alt'] < fdr_tr:
    #     return 'purple'
    if row['log_fc'] * row['log_pv'] > 0:
        return blue_color
    else:
        return red_color


if __name__ == '__main__':

    top10_names = [
       'CEBPB_HUMAN.tsv',
       'ESR1_HUMAN.tsv',
       'SNAI2_HUMAN.tsv',
       'CREB1_HUMAN.tsv',
       'ANDR_HUMAN.tsv',
       'FOXA1_HUMAN.tsv',
       'SPI1_HUMAN.tsv',
       'DUX4_HUMAN.tsv',
       'NRF1_HUMAN.tsv',
       'CTCF_HUMAN.tsv', 'CTCF_Mathelier.tsv']

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
    fc_tr = 2
    fdr_tr_mat = 0.0000000000000000023
    fdr_tr = 0.05

    # blue_color = '#1B7837'
    # red_color = '#762A83'

    blue_color = '#005AB5'
    red_color = '#DC3220'

    # blue_color = '#56B4E9'
    # blue_color = '#0072B2'
    # red_color = '#D55E00'
    grey_color = '#CCCCCC'
    bar_width = 0.8
    bar_alpha = 0.5
    point_lw = 0
    point_alpha = 0.7
    point_size = 10

    # # Barplot
    # df = pd.read_table(os.path.expanduser("~/DataForFigures/blue_red_stats.tsv"))
    # # df = pd.read_table(os.path.expanduser("~/DataForFigures/blue_red_stats.tsv"))
    # fig, ax = plt.subplots(figsize=(6, 5*2/3))
    # plt.tight_layout(rect=(0.025, 0, 1, 1))
    # df['sum'] = df["red"] + df["blue"]
    # df['part'] = df["blue"] / df['sum']
    # df = df[(df['sum'] >= 230)].sort_values("sum")
    # print(df['name'])
    # x, blue, red, blue_n, ticks = range(len(df.index)), \
    #                               df["part"].tolist(), (df["red"] / df['sum']).tolist(), df["blue"].tolist(), \
    #                               df["name"].apply(lambda x: x.replace("_HUMAN", ""))
    #
    # ax.barh(x, red, left=blue, color=red_color, height=bar_width, linewidth=0, tick_label=ticks, alpha=bar_alpha)
    # ax.barh(x, blue, color=blue_color, height=bar_width, linewidth=0, alpha=bar_alpha)
    # for i in x:
    #     ax.text(x=0.5, y=i - 0.1, s='{}/{:.0f}'.format(blue_n[i], blue_n[i] / blue[i]), va="center", ha="center",
    #             fontdict={"size": 13})
    # ax.set_xlabel('Concordant ASB SNPs fraction')
    # ax.tick_params(axis="y", length=0)
    # # plt.savefig(os.path.expanduser("~/AC_7/AS_Figure_7_barplot.png"), dpi=300)
    # plt.savefig(os.path.expanduser("~/AC_7/AS_Figure_7_barplot.svg"), dpi=300)
    #
    # plt.close(fig)

    # Scatters
    # for name in ['All_TFs.tsv'] + top10_names + ['CTCF_for_comparison.tsv']:
    for name in os.listdir(os.path.expanduser('~/Releases/TF_P-values/TF_P-values/')):
        mat = False
        if name == 'CTCF_for_comparison.tsv':
            mat = True
        #pt = pd.read_table(os.path.expanduser("~/DataForFigures/scatter/{}".format('CTCF_HUMAN.tsv' if mat else name)))
        pt = pd.read_table(os.path.expanduser("~/Releases/TF_P-values/TF_P-values/{}".format('CTCF_HUMAN.tsv' if mat else name)))

        pt = pt[(pt['motif_log_pref'] >= -np.log10(perf_tr)) | (pt['motif_log_palt'] >= -np.log10(perf_tr))]
        if 'Mathelier' not in name:
            pt = pt[~(pt[field + '_alt'].isnull() | pt[field + '_ref'].isnull())]
            pt = pt[~(pt['motif_fc'].isnull())]
            pt['log_pv'] = (np.log10(
                pt[[field + '_ref', field + '_alt']]).min(axis=1)) \
                           * np.sign(pt[field + '_alt'] - pt[field + '_ref'])
            pt['log_fc'] = pt['motif_fc']
        else:
            pt['log_fc'] = pt['motif_fc'] / np.log10(2)  # FIXME

        if mat:
            fdr_tr, fdr_tr_mat = fdr_tr_mat, fdr_tr

        if pt.empty:
            continue
        print(pt.head())
        pt['col'] = pt.apply(lambda x: get_color(x), axis=1)

        blue = len(pt[pt['col'] == blue_color].index)
        red = len(pt[pt['col'] == red_color].index)
        grey = len(pt[pt['col'] == grey_color].index)

        if blue + red == 0:
            continue

        pt_grey = pt[pt['col'] == grey_color]
        pt_red = pt[pt['col'] == red_color]
        pt_blue = pt[pt['col'] == blue_color]

        fig, ax = plt.subplots()
        plt.tight_layout(pad=2.5)
        ax.scatter(x=pt_grey['log_pv'], y=pt_grey['log_fc'], c=grey_color, s=point_size, alpha=point_alpha, lw=point_lw, zorder=2, label=None)
        ax.scatter(x=pt_red['log_pv'], y=pt_red['log_fc'], c=red_color, s=point_size, alpha=point_alpha, lw=point_lw, zorder=3, label='Discordant:\n{} ({:.1f})%'.format(red, red / (blue + red) * 100))
        ax.scatter(x=pt_blue['log_pv'], y=pt_blue['log_fc'], c=blue_color, s=point_size, alpha=point_alpha, lw=point_lw, zorder=4, label='Concordant:\n{} ({:.1f})%'.format(blue, blue / (blue + red) * 100))
        ax.set_xlabel('ASB significance')
        ax.set_ylabel('Motif fold change (Alt vs Ref)')

        xlim = max(np.abs(pt['log_pv'])) * 1.05
        ylim = max(np.abs(pt['log_fc'])) * 1.06

        ax.set_xlim(-xlim, xlim)
        ax.set_ylim(-ylim, ylim)

        ax.set_xticklabels(['{:.0f}'.format(abs(x)) for x in ax.get_xticks()])
        ax.set_yticklabels(['{}'.format(x) for x in ax.get_yticks()])
        ax.tick_params(axis='both', zorder=500)
        # print(xlim, ylim, ax.get_xticks(), [x for x in ax.get_xticklabels()])
        # label = 'blue/red: {}/{}({:.1f}%),\ngrey/all={:.1f}%'.format(blue, red,
        #                                                              100 * blue / (blue + red),
        #                                                              100 * grey / (grey + blue + red))
        label = 'concordant/\ndiscordant:\n{}/{}({:.1f}%)'.format(blue, red, 100 * blue / (blue + red))
        # ax.text(x=-xlim * 0.9, y=0, ha='left', va='center', s='Ref', style='italic', color='#505050')
        # ax.text(x=xlim * 0.9, y=0, ha='right', va='center', s='Alt', style='italic', color='#505050')
        ax.text(x=0.22, y=0.066, ha='left', va='bottom', s='Ref', style='italic', color='0.15',
                transform=plt.gcf().transFigure)
        ax.text(x=0.81, y=0.066, ha='right', va='bottom', s='Alt', style='italic', color='0.15',
                transform=plt.gcf().transFigure)
        # plt.text(x=max(pt['log_pv']) / 6, y=min(pt['log_fc']) / 2, s=label)  #, fontdict={'size': 10})
        # plt.text(x=min(pt['log_pv']) * 4 / 5, y=min(pt['log_fc']) / 2,
        #          s=' P-value treshold: {},\nFC treshold: {}'.format(round(fdr_tr, 2), round(fc_tr, 1)))
        ax.grid(True, zorder=-10)

        # box = ax.get_position()
        # ax.set_position([box.x0, box.y0 + box.height * 0.2, box.width, box.height * 0.8])

        ax.legend(
                  fontsize='small',
                  framealpha=0.5,
                  borderpad=0.3,
                  handletextpad=0,
                  fancybox=False,
                  )

        plt.savefig(os.path.expanduser("~/AC_7/AS_Figure_7_{}_p_tr={:.2f}_fc_tr={:.2f}_fdr_tr={:.2f}_{}.svg".format(
            name.replace('_fc.tsv', ''), perf_tr, fc_tr, fdr_tr, field)), dpi=300)
        # plt.savefig(os.path.expanduser("~/AC_7/AS_Figure_7_{}_p_tr={:.2f}_fc_tr={:.2f}_fdr_tr={:.2f}_{}.png".format(
        #     name.replace('_fc.tsv', ''), perf_tr, fc_tr, fdr_tr, field)), dpi=300)
        plt.close(fig)

        if mat:
            fdr_tr, fdr_tr_mat = fdr_tr_mat, fdr_tr
