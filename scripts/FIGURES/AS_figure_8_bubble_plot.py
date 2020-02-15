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
plt.rcParams['figure.figsize'] = 6, 5
plt.rcParams["legend.framealpha"] = 1
# plt.rcParams['axes.xmargin'] = 0
# plt.rcParams['axes.ymargin'] = 0
plt.rcParams["legend.framealpha"] = 1
# rc('text', usetex=True)

if __name__ == '__main__':
    file_name = 'all_lines_union.tsv'
    cl = 'all_lines no filter'
    cosmics = [1, 4/3, 3/2, 5/3, 2, 5/2, 3, 7/2, 11/3, 4, 9/2, 5, 11/2, 6, 7]
    df_counts = pd.read_table(os.path.expanduser('~/Documents/ASB/Correlation/counts.tsv'))
    conv = {1.3333333333333337: 4/3,
            1.6666666666666663: 5/3,
            2.333333333333333: 7/3,
            2.6666666666666665: 8/3,
            3.333333333333333: 10/3,
            3.666666666666667: 11/3}
    for val in conv:
        df_counts = df_counts.replace(val, conv[val])
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

    # Decorations
    plt.savefig(os.path.expanduser('~/AC_8/AS_Figure_8_Counts_plot_{}.png'.format(cl)), dpi=300)
    plt.savefig(os.path.expanduser('~/AC_8/AS_Figure_8_Counts_plot_{}.svg'.format(cl)), dpi=300)
