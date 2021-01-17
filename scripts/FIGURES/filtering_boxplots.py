import itertools
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.ticker import AutoMinorLocator

sns.set(font_scale=1.55, style="ticks", font="lato", palette=('#505050', '#56B4E9', '#009E73', '#F0E442', '#0072B2',
                                                             '#CC79A7', '#D55E00', '#E69F00'))
# sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})
plt.rcParams['font.weight'] = "medium"
plt.rcParams['axes.labelweight'] = 'medium'
plt.rcParams['figure.titleweight'] = 'medium'
plt.rcParams['axes.titleweight'] = 'medium'
plt.rcParams['figure.figsize'] = 16, 7
plt.rcParams["legend.framealpha"] = 1
plt.rcParams['axes.xmargin'] = 0
plt.rcParams['axes.ymargin'] = 0


if __name__ == '__main__':
    import seaborn as sns
    import matplotlib.pyplot as plt
    df = pd.read_table(os.path.expanduser('~/Desktop/table.txt'))
    for column in df.columns:
        if column == 'BaalChIP' or column == 'median BAD ASB' or column == 'Name':
            continue
        print(column)
        df[column] = df[column] / df['BaalChIP']
    df.rename(columns={
        'Called': 'Total set of SNP calls',
        'Filtered': 'SNPs passing basic coverage filters',
        'Filtered (in BaalChIP exps)': 'SNPs passing basic coverage filters\nin BaalChIP-overlapping data',
        'Called (in BaalChIP exps)': 'SNP calls\nin BaalChIP-overlapping data',
        'Candidates': 'SNPs passing complete ADASTRA filters\nfor candidate ASB sites',
        'FDR 0.25': 'ASBs passing FDR ≤ 0.25',
        'FDR 0.2': 'ASBs passing FDR ≤ 0.2',
        'FDR 0.15': 'ASBs passing FDR ≤ 0.15',
        'FDR 0.1': 'ASBs passing FDR ≤ 0.1',
        'FDR 0.05': 'ASBs passing FDR ≤ 0.05',
    }, inplace=True)
    df.set_index('Name', drop=True, inplace=True)
    order = ['IMR90', 'MCF10', 'GM12892', 'GM12891', 'HL60', 'H1hESC', 'HeLa', 'SKNSH', 'T47D', 'A549', 'GM12878', 'HepG2', 'MCF7', 'K562']
    df = df.reindex(order)

    df = df.transpose().head(11)
    df['melt_col'] = df.index
    m_df = df.melt(id_vars='melt_col')
    m_df = m_df[m_df['melt_col'] != 'Name']
    print(m_df.columns)
    mks = itertools.cycle(['x', 'o', '+', '^'])
    # markers = [next(mks) for i in m_df["Name"].unique()]
    markers = ['s', 'o', 'x', 'x', 'o', '^', 's', '*', 'x', 'o', 's', '^', '*', '+']
    print(markers)
    hs = itertools.cycle(['#E69F00', '#009E73', '#CC79A7', '#0072B2', '#D55E00'])
    # hues = [next(hs) for i in m_df["Name"].unique()]
    hues = ['#0072B2'] * 3 + ['#E69F00'] * 5 + ['#009E73'] * 6
    hue_dict = dict(zip(['IMR90', 'MCF10', 'GM12892', 'GM12891', 'HL60', 'H1hESC', 'HeLa', 'SKNSH', 'T47D', 'A549', 'GM12878', 'HepG2', 'MCF7', 'K562'], hues))
    m_df['hues'] = m_df.apply(lambda row: hue_dict[row['Name']], axis=1)
    ax1 = sns.pointplot(y='melt_col', x='value', hue='Name', col='hues', data=m_df, markers=markers, orient='h', size=5,
                        palette=hues,
                      order=['Total set of SNP calls',
                             # 'SNP calls\nin BaalChIP-overlapping data',
                             'SNPs passing basic coverage filters',
                             # 'Passed basic coverage filters\nin BaalChIP-overlapping data',
                             'SNPs passing complete ADASTRA filters\nfor candidate ASB sites',
                             'ASBs passing FDR ≤ 0.25',
                             'ASBs passing FDR ≤ 0.15',
                             'ASBs passing FDR ≤ 0.05'],
                        linestyles='dotted',
                        scale=0.8,
                      )
    plt.tight_layout()
    plt.grid(which='both')
    # ax2 = sns.boxplot(x='melt_col', y='value', data=m_df,
    #                   order=['BaalChIP',
    #                          'Called (in BaalChIP exps)',
    #                          'Filtered (in BaalChIP exps)'])
    plt.ylabel(None)
    plt.xlabel('Fraction of SNPs from the BaalChIP ASB set at the stages of ADASTRA pipeline')
    ax1.xaxis.set_minor_locator(AutoMinorLocator(n=2))
    plt.tick_params(
        axis='x',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom=False,  # ticks along the bottom edge are off
        top=False)  # ticks along the top edge are off
    plt.tick_params(
        axis='y',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        left=False,  # ticks along the bottom edge are off
        right=False)  # ticks along the top edge are off

    # Shrink current axis by 20%
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    legend = ax1.legend(
        loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, fancybox=False,
        fontsize='small',
        borderpad=1.2,
        title='# of SNPs'
    )
    # for i in range(14):
    #     legend.get_texts()[i].set_text(name_dict[i])
    legend.get_frame().set_edgecolor('black')

    plt.savefig(os.path.expanduser('~/Desktop/boxplot.png'), dpi=300)
    plt.show()


