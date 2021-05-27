import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, ticker
import seaborn as sns

sns.set(font_scale=1.4, style="ticks", font="lato",
        # palette=('#56B4E9', '#009E73', '#F0E442'))
        # palette=('#7570b3', '#d95f02', '#1b9e77'))
        palette=('#E69F00', '#009E73', '#D55E00', '#F0E442', '#CC79A7', '#56B4E9'))
# palette=('#1f77b4', '#2ca02c', '#ff7f0e'))
sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})
plt.rcParams['font.weight'] = "medium"
plt.rcParams['axes.labelweight'] = 'medium'
plt.rcParams['figure.titleweight'] = 'medium'
plt.rcParams['axes.titleweight'] = 'medium'
plt.rcParams['figure.figsize'] = 6, 5
plt.rcParams["legend.framealpha"] = 0.5
plt.rcParams['axes.xmargin'] = 0
plt.rcParams['axes.ymargin'] = 0
plt.rcParams["legend.framealpha"] = 1

cells = ('K562', 'MCF7', 'A549', 'HCT116', '22RV1', 'Other')
names = ('K562__myelogenous_leukemia_',
         "MCF7__Invasive_ductal_breast_carcinoma_",
         'A549__lung_carcinoma_',
         'HCT-116__colon_carcinoma_',
         '22RV1__prostate_carcinoma_',
         '')

state_ss = ['int_6', 'full_5_and_6', 'full_6_but_1.33', 'full_6_but_2.5', 'full_6']


def get_cell(row):
    for cell, name in zip(cells, names):
        if row['#cell_line'] == name:
            return cell
    return 'Other'


df = pd.read_table(os.path.expanduser("D:\Sashok/Desktop/susan_BAD/cor_stats_test.tsv"))
df['cell'] = df.apply(get_cell, axis=1)
for state_s in state_ss:
    for i in range(3, 6):
        model = "cor_by_snp_CAIC@{}@{}".format(state_s, i)
        df = df.dropna(subset=[model])
x = sorted(list(df['total_snps'].unique()))

Y = {}
for cell in cells:
    cell_df = df[df['cell'] == cell]
    counts = [cell_df[cell_df['total_snps'] == count]['total_snps'] if not cell_df[cell_df['total_snps'] == count].empty else 0 for count in x ]
    counts = [sum(c) if not isinstance(c, int) else c for c in counts]
    Y[cell] = np.cumsum(counts)

fig, ax = plt.subplots()
fig.tight_layout(pad=2)
for cell in cells:
    ax.plot(x, Y[cell], label=cell)
ax.axvline(x=0, color='#505050', linestyle='--')
ax.axhline(y=0, color='#505050', linestyle='--')
ax.legend(loc='upper left', handletextpad=0.3, handlelength=1)
ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlim(min(x), max(x))
ax.set_ylim(min(x), max(Y['K562']))
ax.axvline(x=2000, color='#50505050', ls='--')
ax.axvline(x=60000, color='#50505050', ls='--')
ax.grid(True)

ax.set_ylabel("Cumulative number of SNPs")
ax.set_xlabel("Number of SNPs in a group of datasets")

plt.title('SNP count distribution')

plt.savefig(os.path.expanduser('D:\Sashok/Desktop/susan_BAD/snp_dist.png'), dpi=300)
plt.close(fig)
