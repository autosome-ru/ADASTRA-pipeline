# import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, ticker
import seaborn as sns
from scipy import stats as st


def make_binom_density(n, p):
    binom = np.zeros(n + 1)
    if p != 0.5:
        f1 = st.binom(n, p).pmf
        f2 = st.binom(n, 1 - p).pmf
        binom_norm = 1 - sum(0.5 * (f1(k) + f2(k)) for k in [0, 1, 2, 3, 4, n - 4, n - 3, n - 2, n - 1, n])
        for k in range(5, n - 4):
            binom[k] = 0.5 * (f1(k) + f2(k)) / binom_norm
    else:
        f = st.binom(n, p).pmf
        binom_norm = 1 - sum(f(k) for k in [0, 1, 2, 3, 4, n - 4, n - 3, n - 2, n - 1, n])
        for k in range(5, n - 4):
            binom[k] = f(k) / binom_norm
    return binom


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



lw = 1.25
states = [1, 2, 3, 4, 5, 3/2]

for max_c, cells, covs, BAD, p, suffix in [
    # [150, 'K562', [30, 60], 0, 1/3, '_big'],
    # [100, 'diploid', [15, 30], 0, 1/2, '_small'],
    # [50, 'HCT116', [15, 30], 0, 1/2, '_small'],
] + [[150, 'all', [30, 60], B, 0, '_big'] for B in states] + [[50, 'all', [15, 30], B, 0, '_small'] for B in states + [0]] + \
        [[50, 'all', [15, 30], 0, 1/2, '_small'], [150, 'all', [30, 60], 0, 1/2, '_big']]:
    if BAD != 0:
        p = 1/(BAD + 1)
    t = pd.read_table(os.path.expanduser(
        '~/DataForFigures/bias_stats_BAD{:.1f}{}.tsv'.format(BAD if BAD else 0,
                                                              {'all': '',
                                                               'K562': '_k562',
                                                               'diploid': '_esc',
                                                               'HCT116': '_hct116'}[
                                                                  cells])))

    t.columns = ['alt', 'ref', 'count']

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

    t.columns = ['Alternative allele read count', 'Reference allele read count', 'count']
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

    ax1.hlines(y=max_c + 1 - 5, xmin=0, xmax=max_c + 1 - 5, colors=['black', ], linewidth=lw*2)
    ax1.vlines(x=0, ymin=0, ymax=max_c + 1 - 5, colors=['black', ], linewidth=lw*2)
    ax1.hlines(y=0, xmin=0, xmax=max_c + 1 - 5, colors=['black', ], linewidth=lw*2)
    ax1.vlines(x=max_c + 1 - 5, ymin=0, ymax=max_c + 1 - 5, colors=['black', ], linewidth=lw*2)

    for cov in covs:
        ax1.plot([0, cov], [max_c - cov + 0.5, max_c + 0.5], linestyle='dashed', linewidth=lw, color='black')


    for cov, ax in zip(covs, (ax2, ax3)):
        counts_array = np.array(arrays[cov])
        total_snps = counts_array.sum()

        x = list(range(cov + 1))
        sns.barplot(x=x, y=counts_array / total_snps, ax=ax, color='C1')

        ax.xaxis.set_major_locator(ticker.FixedLocator(np.arange(0, len(x), 5)))
        ax.xaxis.set_major_formatter(ticker.FixedFormatter(x[::5]))
        ax.tick_params(axis="x", rotation=0)

        current_density = list(make_binom_density(cov, p=p))

        slope, intercept, r_value, p_value, std_err = st.linregress(x[5: -5], counts_array[5: -5] / total_snps)
        print(slope, intercept)
        ax.plot(x, np.array(x) * slope + intercept, color='#DC3220')

        ax.plot(sorted(x + [5, cov - 5]), [0] + current_density + [0], color='#4d004b')
        ax.set_ylim(0, max(max(current_density), max(counts_array / total_snps)) * 1.05)
        ax.set_xlabel('Reference allele read count')

    # plt.savefig(os.path.expanduser('~/AC_3/Figure_AS_3_{}.svg'.format(cells)))
    # plt.savefig(os.path.expanduser('~/AC_3/Figure_AS_3_{}.png'.format(cells)))
    plt.savefig(os.path.expanduser('~/AC_3/Figure_AS_3_{}_{:.2f}{}.svg'.format(cells, BAD, suffix)))
    plt.savefig(os.path.expanduser('~/AC_3/Figure_AS_3_{}_{:.2f}{}.png'.format(cells, BAD, suffix)))
    # plt.show()
