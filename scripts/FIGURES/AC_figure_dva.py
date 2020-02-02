import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, ticker
import seaborn as sns

sns.set(font_scale=1.4, style="ticks", font="lato", palette=('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2',
                                                             '#D55E00', '#CC79A7'))
# sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})
plt.rcParams['font.weight'] = "medium"
plt.rcParams['axes.labelweight'] = 'medium'
plt.rcParams['figure.titleweight'] = 'medium'
plt.rcParams['axes.titleweight'] = 'medium'
plt.rcParams['figure.figsize'] = 6, 5
plt.rcParams["legend.framealpha"] = 1
plt.rcParams['axes.xmargin'] = 0
plt.rcParams['axes.ymargin'] = 0
plt.rcParams["legend.framealpha"] = 1

max_c = 150
min_c = 0
lw = 1.25

t = pd.read_table(os.path.expanduser('~/all_snps_statistics.tsv'))
t = t[(t['ref'] >= min_c) & (t['alt'] >= min_c)]
t = t[(t['ref'] <= max_c) & (t['alt'] <= max_c)]
t['ref'] = t['ref'].astype(int)
t['alt'] = t['alt'].astype(int)
print(t['count'].sum(axis=0))
t.columns = ['Reference allele read count', 'Alternative allele read count', 'count']
t = t.pivot('Alternative allele read count', 'Reference allele read count', 'count')
t.sort_index(ascending=False, inplace=True)
t.fillna(0, inplace=True)

for l in range(min_c, max_c + 1):
    for k in range(min_c, l + 1):
        if t[k][l] + t[l][k] == 0:
            t[k][l] = 0
            t[l][k] = 0
        else:
            t[k][l], t[l][k] = ((t[k][l] - t[l][k]) / (t[k][l] + t[l][k]) * np.log10(t[k][l] + t[l][k] + 1),
                                (t[l][k] - t[k][l]) / (t[k][l] + t[l][k]) * np.log10(t[k][l] + t[l][k] + 1))

fig, ax = plt.subplots()
fig.tight_layout(pad=2)
sns.heatmap(t, cmap="PRGn", vmin=-5, vmax=5)

cbar = ax.collections[0].colorbar
cbar.set_ticks(np.arange(-5, 6, 1))

# SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
# SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
#
# cbar.set_ticklabels(["-10⁵", "-10⁴", "-10³", "-10²", "-10¹", "0", "10¹", "10²", "10³", "10⁴", "10⁵"])

if max_c <= 50:
    div = 5
else:
    div = 10

ax.yaxis.set_major_locator(ticker.FixedLocator(np.arange(0, max_c - 5 + 1, div) + 0.5))
ax.yaxis.set_major_formatter(ticker.FixedFormatter(range(0, max_c + 1)[::-div]))
ax.tick_params(axis="y", rotation=0)

ax.xaxis.set_major_locator(ticker.FixedLocator(np.arange(0, max_c + 1, div) + 0.5))
ax.xaxis.set_major_formatter(ticker.FixedFormatter(range(0, max_c + 1)[::div]))
ax.tick_params(axis="x", rotation=-90)

plt.hlines(y=max_c + 1, xmin=0, xmax=max_c + 1, colors=['black', ], linewidth=lw*2 + 0.5)
plt.vlines(x=0, ymin=0, ymax=max_c + 1, colors=['black', ], linewidth=lw*2)


plt.hlines(y=0, xmin=5, xmax=max_c + 1, colors=['#AAAAAA', ], linewidth=lw*2)
plt.vlines(x=max_c + 1, ymin=0, ymax=max_c + 1 - 5, colors=['#AAAAAA', ], linewidth=lw*2 + 0.5)
plt.hlines(y=max_c + 1 - 5, xmin=5-lw/10, xmax=max_c + 1, colors=['#AAAAAA', ], linewidth=lw)
plt.vlines(x=5, ymin=0, ymax=max_c + 1 - 5 +lw/10, colors=['#AAAAAA', ], linewidth=lw)
plt.savefig(os.path.expanduser('~/AC2/Figure_AS_2_290120.svg'), dpi=300)
plt.show()
