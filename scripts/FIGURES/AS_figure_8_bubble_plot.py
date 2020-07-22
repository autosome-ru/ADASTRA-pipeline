import sys
import os
import numpy as np
import pandas as pd
from matplotlib import rc

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.helpers import states
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns

sns.set(font_scale=1, style="ticks", font="lato")
sns.set_palette((
'#56B4E9',
'#0072B2',
'#009E73',
'#FFFFFF',
'#E69F00',
'#F0E442',
'#D55E00',
'#FFFFFF',
'#FFFFFF',
'#999999',
'#FFFFFF',
'#505050',
'#FFFFFF',
'#CC79A7',
'#FFFFFF'))
sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})
plt.rcParams['font.weight'] = "medium"
plt.rcParams['axes.labelweight'] = 'medium'
plt.rcParams['figure.titleweight'] = 'medium'
plt.rcParams['axes.titleweight'] = 'medium'
plt.rcParams['figure.figsize'] = 6.5, 5
plt.rcParams["legend.framealpha"] = 1
# plt.rcParams['axes.xmargin'] = 0
# plt.rcParams['axes.ymargin'] = 0
plt.rcParams["legend.framealpha"] = 1
# rc('text', usetex=True)

if __name__ == '__main__':
    file_name = 'all_lines_union.tsv'
    cl = 'all_lines no filter'
    cosmics = [1, 4/3, 3/2, 5/3, 2, 5/2, 3, 7/2, 11/3, 4, 9/2, 5, 11/2, 6, 7]
    df_counts = pd.read_table(os.path.expanduser('~/DataForFigures/counts.tsv'))
    conv = {1.3333333333333337: 4/3,
            1.6666666666666663: 5/3,
            2.333333333333333: 7/3,
            2.6666666666666665: 8/3,
            3.333333333333333: 10/3,
            3.666666666666667: 11/3}
    for val in conv:
        df_counts = df_counts.replace(val, conv[val])
    for f in states:
        print('{:.2f}: {:.2f}'.format(f, df_counts[df_counts['COSMIC'].isin([f])]['counts'].sum() / df_counts['counts'].sum() * 100))
    print('{}: {:.2f}'.format('other', df_counts[~df_counts['COSMIC'].isin(states)]['counts'].sum() / df_counts['counts'].sum() * 100))
    df_counts = df_counts[df_counts['COSMIC'].isin(cosmics)]
    # Draw Stripplot
    fig, ax = plt.subplots()
    plt.tight_layout(pad=2.5)
    ax.set_yticks(list(set(df_counts.BAD)), minor=False)
    ax.set_xticks(list(set(df_counts.COSMIC)), minor=False)
    ax.yaxis.grid(True, which='major')
    ax.xaxis.grid(True, which='major')
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    print(sorted(list(set(df_counts['COSMIC']))))
    ax.set_xticks([1, 4/3, 3/2, 5/3, 2, 5/2, 3, 7/2, 11/3, 4, 9/2, 5, 11/2, 6, 7])
    ax.set_xticklabels(
        ['1', '4/3', '3/2', '5/3', '2', '5/2', '3', '7/2', '11/3', '4', '9/2', '5', '11/2', '6', '>6']
    )
    ax.set_yticklabels(
        ['1', '4/3', '3/2', '2', '5/2', '3', '4', '5', '6']
    )
    ax.set_yticks(states)
    ax.tick_params(axis="y", rotation=0)
    ax.tick_params(axis="x", rotation=-90)

    for i, c_bad in enumerate(cosmics):
        if i == 0:
            leg = 'brief'
        else:
            leg = False
        sns.scatterplot(df_counts[df_counts['COSMIC'] == c_bad]['COSMIC'],
                        df_counts[df_counts['COSMIC'] == c_bad]['BAD'],
                        size=df_counts.counts ** 0.5, sizes=(2, 800),
                        alpha=0.7,
                        linewidth=0.5,
                        # color='C' + str(i),
                        edgecolor='#505050',
                        legend=False)

    ax.set_axisbelow(True)
    ax.set_ylabel('Segmentation BAD')
    ax.set_xlabel('COSMIC BAD')

    print(max(df_counts.counts ** 0.5), min(df_counts.counts ** 0.5))

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    legend_elements = [
        plt.scatter([-1], [-1], s=(10**(5/2)-1) / 4897.98 * 798 + 2,
                        alpha=0.7,
                        linewidth=0.5,
                        color='#505050',
                        edgecolor='#505050',
                        label='1 x 10⁵'),
        plt.scatter([-1], [-1], s=(10**(6/2)-1) / 4897.98 * 798 + 2,
                        alpha=0.7,
                        linewidth=0.5,
                        color='#505050',
                        edgecolor='#505050',
                        label='1 x 10⁶'),
        plt.scatter([-1], [-1], s=(10**(7/2)-1) / 4897.98 * 798 + 2,
                        alpha=0.7,
                        linewidth=0.5,
                        color='#505050',
                        edgecolor='#505050',
                        label='1 x 10⁷'),
    ]

    # Put a legend to the right of the current axis
    legend = ax.legend(
        loc='center left', bbox_to_anchor=(1, 0.5), frameon=False, fancybox=False,
        handles=legend_elements,
        fontsize='small',
        labelspacing=2.5,
        borderpad=1.2,
        title='# of SNPs'
    )
    legend.get_frame().set_edgecolor('black')

    ax.set_xlim(0.67, 6.3)
    ax.set_ylim(0.67, 6.25)

    plt.savefig(os.path.expanduser('~/AC_8/AS_Figure_8_Counts_plot_{}.png'.format(cl)), dpi=300)
    plt.savefig(os.path.expanduser('~/AC_8/AS_Figure_8_Counts_plot_{}.svg'.format(cl)), dpi=300)
