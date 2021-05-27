import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, ticker
import seaborn as sns

sns.set(font_scale=1.4, style="ticks", font="lato",
        # palette=('#56B4E9', '#009E73', '#F0E442'))
        # palette=('#7570b3', '#d95f02', '#1b9e77'))
        palette=('#56B4E9', '#E69F00', '#009E73', '#D55E00', '#CC79A7', '#F0E442'))
# palette=('#1f77b4', '#2ca02c', '#ff7f0e'))
sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})
plt.rcParams['font.weight'] = "medium"
plt.rcParams['axes.labelweight'] = 'medium'
plt.rcParams['figure.titleweight'] = 'medium'
plt.rcParams['axes.titleweight'] = 'medium'
plt.rcParams['figure.figsize'] = 16, 8
plt.rcParams["legend.framealpha"] = 0.5
plt.rcParams['axes.xmargin'] = 0
plt.rcParams['axes.ymargin'] = 0
plt.rcParams["legend.framealpha"] = 1


df_stats = pd.read_table(os.path.expanduser("D:\Sashok/Desktop/susan_BAD/stats.tsv"))
cells = ('K562', 'MCF7', 'A549', 'HCT116', '22RV1', 'Other')
markers = ('s', '*', 'v', '1', 'o')
state_ss = ['int_6', 'full_5_and_6', 'full_6_but_1.33', 'full_6_but_2.5', 'full_6']
colors = ('#E69F00', '#009E73', '#D55E00', '#CC79A7', '#F0E442', '#56B4E9')
fig, ax = plt.subplots()
for i, (cell_sign, color) in enumerate(zip(cells, colors), 1):
    for mult in range(3, 6):
        for state_s, marker in zip(state_ss, markers):
            for j, size in zip(range(3), ('small', 'med', 'big')):
                row = df_stats[(df_stats['States'] == state_s) & (df_stats['Multiplier'] == mult)]
                ax.scatter([i + (mult-4)/9 + (j-1)/3], row['{}_{}_Q3'.format(cell_sign, size)],
                           s=40 * (mult - 2), zorder=10, alpha=0.5, marker=marker, color=color)

for i in range(1, len(cells)):
    ax.axvline(x=i+1/2, color='#505050')
for i in range(1, len(cells) + 1):
    for j, text, hint in zip(range(3), ('small', 'med', 'big'), ('<2K', '<60K', '≥60K')):
        if i == 1:
            ax.text(y=0.025, x=i + (j - 1) / 3, s=hint+'\n'+text, ha='center', va='bottom', color='#50505050')
        else:
            ax.text(y=0.025, x=i + (j - 1) / 3, s=text, ha='center', va='bottom', color='#50505050')
        if j != 1:
            ax.axvline(x=i + (j-1)/6, color='#50505050', ls='--')
ax.set_ylim(0, 1)
ax.set_xlim(1/2, len(cells) + 1/2)
ax.grid(axis='y')
ax.set_xticks(list(range(1, len(cells) +1)))
ax.set_xticklabels(cells)
ax.set_ylabel("Kendall's tau (Segmentation, COSMIC), Q3")
ax.set_xlabel("Cell line name")
ax.set_title('3-d quartile of correlation with COSMIC\ndepending on states set and CAIC multiplier')

# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
legend_elements = [
                      plt.scatter([-1], [-1], color='black', alpha=0.5, marker=marker) for marker in markers
                  ]
legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='BAD / states', frameon=True, fancybox=False,
                   labels=state_ss, ncol=1, handles=legend_elements, handlelength=1)
legend.get_frame().set_edgecolor('black')

plt.savefig(os.path.expanduser('D:\Sashok/Desktop/susan_BAD/cor_sum_plot.png'), dpi=300)

