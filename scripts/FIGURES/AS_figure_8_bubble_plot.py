import sys
import os
import numpy as np
import pandas as pd

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.FIGURES import style_config
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns

sns.set(font_scale=1.4, style="ticks", font="lato", palette=('#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00',
                                                             '#ffff33', '#a65628', '#f781bf', '#999999'))
sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})
plt.rcParams['font.weight'] = "medium"
plt.rcParams['axes.labelweight'] = 'medium'
plt.rcParams['figure.titleweight'] = 'medium'
plt.rcParams['axes.titleweight'] = 'medium'
plt.rcParams['figure.figsize'] = 8, 5
plt.rcParams["legend.framealpha"] = 1
# plt.rcParams['axes.xmargin'] = 0
# plt.rcParams['axes.ymargin'] = 0
plt.rcParams["legend.framealpha"] = 1

if __name__ == '__main__':
    file_name = 'all_lines_union.tsv'
    cl = 'all_lines no filter'
    df_counts = pd.read_table(os.path.expanduser('~/Documents/ASB/Correlation/counts.tsv'))
    # Draw Stripplot
    fig, ax = plt.subplots()
    plt.tight_layout(pad=1.5)
    ax.set_yticks(list(set(df_counts.BAD)), minor=False)
    ax.set_xticks(list(set(df_counts.COSMIC)), minor=False)
    ax.yaxis.grid(True, which='major')
    ax.xaxis.grid(True, which='major')
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    print(sorted(list(set(df_counts['COSMIC']))))
    ax.set_xticks([1, 5/4, 4/3, 3/2, ])
    ax.tick_params(axis="y", rotation=90)

    sns.scatterplot(df_counts.COSMIC, df_counts.BAD, size=df_counts.counts ** 0.5, sizes=(2, 800),
                    hue=df_counts.COSMIC, palette=sns.color_palette('husl', len(set(df_counts.COSMIC))),
                    legend=False)

    ax.set_axisbelow(True)

    # Decorations
    plt.savefig(os.path.expanduser('~/AC_8/AS_Figure_8_Counts_plot_{}.png'.format(cl)), dpi=300)
