import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, ticker
import seaborn as sns
from scipy import stats as st

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.helpers import states, read_weights


def make_negative_binom_density(r, p, w, size_of_counts):
    negative_binom_density_array = np.zeros(size_of_counts + 1, dtype=np.float128)
    dist1 = st.nbinom(r, p)
    f1 = dist1.pmf
    cdf1 = dist1.cdf
    dist2 = st.nbinom(r, 1 - p)
    f2 = dist2.pmf
    cdf2 = dist2.cdf
    negative_binom_norm = (cdf1(size_of_counts) - cdf1(4)) * w + \
                          (cdf2(size_of_counts) - cdf2(4)) * (1 - w)

    for k in range(5, size_of_counts + 1):
        negative_binom_density_array[k] = (w * f1(k) + (1 - w) * f2(k)) / negative_binom_norm
    return negative_binom_density_array


sns.set(font_scale=1.55, style="ticks", font="lato", palette=('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2',
                                                             '#D55E00', '#CC79A7'))
# sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})
plt.rcParams['font.weight'] = "medium"
plt.rcParams['axes.labelweight'] = 'medium'
plt.rcParams['figure.titleweight'] = 'medium'
plt.rcParams['axes.titleweight'] = 'medium'
plt.rcParams['figure.figsize'] = 18, 5
plt.rcParams["legend.framealpha"] = 1
plt.rcParams['axes.xmargin'] = 0
plt.rcParams['axes.ymargin'] = 0
plt.rcParams["legend.framealpha"] = 1

max_c = 50
lw = 1.25
cells = 'LoVo'
p = 1 / 2
covs = [15, 30]

t = pd.read_table(os.path.expanduser('~/{}_snps_statistics.tsv'.format(cells)))
t = t[(t['ref'] >= 5) & (t['alt'] >= 5)]
t = t[(t['ref'] <= max_c) & (t['alt'] <= max_c)]
for count in range(5, max_c + 1):
    if not t[(t['ref'] == count) & (t['alt'] == 5)]['count'].tolist():
        t = t.append(pd.DataFrame({'ref': [count], 'alt': [5], 'count': [0]}))
    if not t[(t['ref'] == 5) & (t['alt'] == count)]['count'].tolist():
        t = t.append(pd.DataFrame({'ref': [5], 'alt': [count], 'count': [0]}))
t.reset_index(inplace=True, drop=True)
t['ref'] = t['ref'].astype(int)
t['alt'] = t['alt'].astype(int)
t['count'] = t['count'].astype(int)
print(t['count'].sum(axis=0))

arrays = {}
for cov in covs:
    t_slice = t[t['ref'] + t['alt'] == cov][['ref', 'count']]
    arrays[cov] = [t_slice[t_slice['ref'] == k]['count'].tolist()[0]
                   if t_slice[t_slice['ref'] == k]['count'].tolist()
                   else 0
                   for k in range(cov + 1)]

t.columns = ['Reference allele read count', 'Alternative allele read count', 'count']
t = t.pivot('Alternative allele read count', 'Reference allele read count', 'count')
t.sort_index(ascending=False, inplace=True)
t.fillna(0, inplace=True)

t = np.log10(t + 1)

fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
fig.tight_layout(pad=1.5)
sns.heatmap(t, cmap="BuPu", ax=ax1, vmin=0, vmax=4)

cbar = ax1.collections[0].colorbar
cbar.set_ticks(np.arange(0, 5, 1))
cbar.set_ticklabels(["10⁰", "10¹", "10²", "10³", "10⁴"])

if max_c <= 50:
    div = 5
elif max_c <= 100:
    div = 10
elif max_c % 30 == 0:
    div = 30
else:
    div = 50

ax1.yaxis.set_major_locator(ticker.FixedLocator(np.arange(0, max_c - 5 + 1, div) + 0.5))
ax1.yaxis.set_major_formatter(ticker.FixedFormatter(range(div - 5, max_c + 1)[::-div]))
ax1.tick_params(axis="y", rotation=0)

ax1.xaxis.set_major_locator(ticker.FixedLocator(np.arange(div - 5, max_c + 1, div) + 0.5))
ax1.xaxis.set_major_formatter(ticker.FixedFormatter(range(div, max_c + 1)[::div]))
ax1.tick_params(axis="x", rotation=0)

ax1.hlines(y=max_c + 1 - 5, xmin=0, xmax=max_c + 1 - 5, colors=['black', ], linewidth=lw*2 + 1)
ax1.vlines(x=0, ymin=0, ymax=max_c + 1 - 5, colors=['black', ], linewidth=lw*2)
ax1.hlines(y=0, xmin=0, xmax=max_c + 1 - 5, colors=['black', ], linewidth=lw*2)
ax1.vlines(x=max_c + 1 - 5, ymin=0, ymax=max_c + 1 - 5, colors=['black', ], linewidth=lw*2 + 1)

for cov in covs:
    ax1.plot([0, cov], [max_c - cov + 0.5, max_c + 0.5], linestyle='dashed', linewidth=lw, color='black')

# Params
r_dict, w_dict, gof_dict = read_weights()
fix_c = 20
main_allele = "alt"
fixed_allele = "ref" if main_allele == "alt" else "alt"
BAD = 1
if BAD not in states:
    print("Sanya hueviy BAD")
stats = pd.read_table(os.path.expanduser('~/fixed_alt_bias_statistics_BAD={:.1f}.tsv'.format(BAD)))
stats_filtered = stats[stats['{}_counts'.format(fixed_allele)] == fix_c]
max_cover_in_stats = max(stats_filtered['{}_counts'.format(main_allele)])
counts_array = np.zeros(max_cover_in_stats + 1, dtype=np.int64)
for index, row in stats_filtered.iterrows():
    k, SNP_counts = row['{}_counts'.format(main_allele)], row['counts']
    counts_array[k] = SNP_counts

total_snps = counts_array[0:max_cover_in_stats + 1].sum()
fig, ax = plt.subplots(figsize=(10, 8))
x = list(range(max_cover_in_stats + 1))
sns.barplot(x=x,
            y=counts_array[0:max_cover_in_stats + 1] / total_snps, ax=ax)
ax.xaxis.set_major_locator(ticker.FixedLocator(np.arange(0, len(x), 5)))
ax.xaxis.set_major_formatter(ticker.FixedFormatter(x[::5]))
ax.tick_params(axis="x", rotation=90)

r, w, gof = (r_dict[fixed_allele][BAD][fix_c],
             w_dict[fixed_allele][BAD][fix_c],
             gof_dict[fixed_allele][BAD][fix_c])
current_density = make_negative_binom_density(r, BAD, w, fix_c)

asb = counts_array[5:].sum()

label = 'negative binom fit for {}r_dict, w_dict, ' \
        'gof_dict = read_weights()\ntotal observations: {}\nr={:.2f}, p={:.2f}, q={}, w={:.2f}\ngof={:.4}'.format(
    main_allele, total_snps, r, BAD, w, asb, gof)
plt.plot(list(range(fix_c + 1)), current_density)
plt.text(s=label, x=0.65 * fix_c, y=max(current_density) * 0.6)
plt.title('scaled ref: fixed_{}={}, BAD={:.1f}, 2 params'.format(fixed_allele, fix_c, BAD))
plt.savefig(os.path.expanduser(
    '~/fixed_alt/abcd/scaled_2params_q15-q95_{}_BAD={:.1f}_fixed_{}.png'.format(fixed_allele, BAD, fix_c)))
plt.close(fig)
