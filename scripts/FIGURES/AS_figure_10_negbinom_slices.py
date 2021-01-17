import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, ticker
import seaborn as sns
from scipy import stats as st

states = [1, 1.5, 2, 3, 4, 5]

BAD_dict = {1: '1', 2: '2', 3: '3', 4: '4', 5: '5', 6: '6', 4/3: '4/3', 3/2: '3/2', 5/2: '5/2'}


def make_negative_binom_density(r, p, w, size_of_counts):
    negative_binom_density_array = np.zeros(size_of_counts + 1, dtype=np.float_)
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


def local_read_weights():
    r = {}
    w = {}
    gof = {}
    for fixed_allele in ('ref', 'alt'):
        r[fixed_allele] = {}
        w[fixed_allele] = {}
        gof[fixed_allele] = {}
        for BAD in states:
            precalc_params_path = os.path.expanduser('~/DataForFigures/NBweights_{}_BAD={:.1f}.npy'.format(fixed_allele, BAD))
            coefs_array = np.load(precalc_params_path)
            r[fixed_allele][BAD] = coefs_array[:, 0]
            w[fixed_allele][BAD] = coefs_array[:, 1]
            gof[fixed_allele][BAD] = coefs_array[:, 3]
            first_bad_gof = min(x for x in range(len(gof[fixed_allele][BAD])) if gof[fixed_allele][BAD][x] > 0.05)
            gof[fixed_allele][BAD][first_bad_gof:] = 1
            r[fixed_allele][BAD][first_bad_gof:] = 0
            w[fixed_allele][BAD][first_bad_gof:] = 1
    return r, w, gof


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
cells = 'all'
covs = [10, 10]
fixs = ['ref', 'alt']
color_maxlog = 4
r_dict, w_dict, gof_dict = local_read_weights()

for BAD in states:
    p = 1 / (BAD + 1)
    if max(covs) * BAD >= 50:
        max_c = 100
    else:
        max_c = 50

    t = pd.read_table(os.path.expanduser(
        '~/DataForFigures/bias_stats_BAD{:.1f}{}.tsv'.format(BAD,
                                                              {'all': '',
                                                               'K562': '_k562',
                                                               'diploid': '_esc'}[
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
    sns.heatmap(t, cmap="BuPu", ax=ax1, vmin=0, vmax=color_maxlog)

    cbar = ax1.collections[0].colorbar
    cbar.set_ticks(np.arange(0, color_maxlog + 1, 1))
    cbar.set_ticklabels(["10⁰", "10¹", "10²", "10³", "10⁴", "10⁵", "10⁶", "10⁷"])

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

    ax1.hlines(y=max_c + 1 - 5, xmin=0, xmax=max_c + 1 - 5, colors=['black', ], linewidth=lw * 2)
    ax1.vlines(x=0, ymin=0, ymax=max_c + 1 - 5, colors=['black', ], linewidth=lw * 2)
    ax1.hlines(y=0, xmin=0, xmax=max_c + 1 - 5, colors=['black', ], linewidth=lw * 2)
    ax1.vlines(x=max_c + 1 - 5, ymin=0, ymax=max_c + 1 - 5, colors=['black', ], linewidth=lw * 2)

    for cov, fix in zip(covs, fixs):
        if fix == 'alt':
            ax1.axhline(xmin=0, y=max_c - cov + 0.5, linestyle='dashed', linewidth=lw, color='black')
        elif fix == 'ref':
            ax1.axvline(x=cov - 5 + 0.5, ymin=0, linestyle='dashed', linewidth=lw, color='black')

    # Params
    for fix_c, fixed_allele, ax in zip(covs, fixs, (ax2, ax3)):
        main_allele = "ref" if fixed_allele == "alt" else "alt"
        stats = pd.read_table(os.path.expanduser(
            '~/DataForFigures/bias_stats_BAD{:.1f}{}.tsv'.format(BAD,
                                                                  {'all': '',
                                                                   'K562': '_k562',
                                                                   'diploid': '_esc'}[
                                                                      cells])))
        stats_filtered = stats[stats['{}_counts'.format(fixed_allele)] == fix_c].astype(int)
        max_cover_in_stats = max(stats_filtered['{}_counts'.format(main_allele)], default=10)
        counts_array = np.zeros(max_cover_in_stats + 1, dtype=np.int64)
        for index, row in stats_filtered.iterrows():
            k, SNP_counts = row['{}_counts'.format(main_allele)], row['counts']
            counts_array[k] = SNP_counts

        chop_counts_array = np.zeros(max_c + 1)
        chop_counts_array[:min(max_c, max_cover_in_stats) + 1] = counts_array[:min(max_c, max_cover_in_stats) + 1]

        total_snps = counts_array[0:max_cover_in_stats + 1].sum()
        x = list(range(max_c + 1))
        sns.barplot(x=x,
                    y=chop_counts_array / total_snps, ax=ax, color='C1')

        r, w, gof = (r_dict[fixed_allele][BAD][fix_c],
                     w_dict[fixed_allele][BAD][fix_c],
                     gof_dict[fixed_allele][BAD][fix_c])
        print(r, w, gof, BAD, fix_c)
        if r == 0:
            col = 'C6'
            r = fix_c
        else:
            col = '#4d004b'

        current_density = np.zeros(max_c + 1)
        current_density[:min(max_c, max_cover_in_stats) + 1] = \
            make_negative_binom_density(r, p, w, max_cover_in_stats)[:min(max_c, max_cover_in_stats) + 1]

        asb = counts_array[5:].sum()

        label = 'negative binom fit for {}' \
                '\ntotal observations: {}\nr={:.2f}, p={:.2f}, q={}, w={:.2f}\ngof={:.4}'.format(
            main_allele, total_snps, r, p, w, asb, gof)
        ax.plot(sorted(x + [5]), [0] + list(current_density), color=col)
        # ax.text(s=label, x=0.65 * fix_c, y=max(current_density) * 0.6)

        ax.xaxis.set_major_locator(ticker.FixedLocator(np.arange(0, max_c + 1, div)))
        ax.xaxis.set_major_formatter(ticker.FixedFormatter(range(0, max_c + 1)[::div]))
        ax.tick_params(axis="x", rotation=0)
        ax.set_xlabel('{} allele read count'.format({'ref': 'Reference', 'alt': 'Alternative'}[main_allele]))

    # plt.suptitle('Figure_AS_10_{}_{:.2f}'.format(cells, BAD))

    plt.savefig(os.path.expanduser('~/AC_10/Figure_AS_10_{}_{:.2f}.svg'.format(cells, BAD)))
    plt.close(fig)

    # gof vs read cov

    df_alt = pd.read_table(os.path.expanduser('~/DataForFigures/NBweights_alt_BAD{:.1f}.tsv'.format(BAD)))
    df_ref = pd.read_table(os.path.expanduser('~/DataForFigures/NBweights_ref_BAD{:.1f}.tsv'.format(BAD)))

    df_ref = df_ref[(df_ref['gof'] > 0) & (df_ref.index <= 50)]
    df_alt = df_alt[(df_alt['gof'] > 0) & (df_alt.index <= 50)]

    fig, ax = plt.subplots(figsize=(6, 5))
    fig.tight_layout(pad=2)

    ax.set_xlim(5, 50)
    ax.set_ylim(0, max(max(df_ref['gof']), max(df_alt['gof'])) * 1.05)
    ax.grid(True)

    ax.axhline(y=0.05, lw=2, linestyle='--', color='#505050')

    ax.scatter(x=df_alt.index, y=df_alt["gof"].tolist(), color='C1', label='Alt')
    ax.scatter(x=df_ref.index, y=df_ref["gof"].tolist(), color='C2', label='Ref')

    ax.set_xlabel('Read count for the fixed allele')
    ax.set_ylabel('Goodness of fit, RMSEA')

    plt.title('BAD={}'.format(BAD_dict[BAD]))

    ax.legend(title='Fixed allele')

    plt.savefig(os.path.expanduser('~/AC_10/Figure_AS_10_gof_{:.2f}.svg'.format(BAD)))
    plt.close(fig)

    # ref bias in r scatter

    df_ref = df_ref[(df_ref['gof'] <= 0.05)]
    df_alt = df_alt[(df_alt['gof'] <= 0.05)]

    fig, ax = plt.subplots(figsize=(6, 5))
    fig.tight_layout(pad=2)

    ax.set_xlim(5, 50)
    y_max = max(max(df_ref['r']), max(df_alt['r']))
    ax.set_ylim(0, y_max * 1.05)
    ax.grid(True)

    ax.plot([5, y_max], [5, y_max], c='grey', label='y=x', linestyle='dashed')
    if BAD == 4 / 3:
        ax.plot([10 * 4 / 3, y_max * 4 / 3], [10, y_max], label='y=3/4 x', c='black', linestyle='dashed')

    ax.scatter(x=df_alt.index, y=df_alt["r"].tolist(), color='C1', label='Alt')
    ax.scatter(x=df_ref.index, y=df_ref["r"].tolist(), color='C2', label='Ref')

    ax.set_xlabel('Read count for the fixed allele')
    ax.set_ylabel('Fitted r value')

    ax.legend(title='Fixed allele')

    plt.title('BAD={}'.format(BAD_dict[BAD]))

    plt.savefig(os.path.expanduser('~/AC_10/Figure_AS_10_r_scatter_{:.2f}.svg'.format(BAD)))
    plt.close(fig)

    # ref bias in r scatter 2

    fig, ax = plt.subplots(figsize=(6, 5))
    fig.tight_layout(pad=2)

    # ax.set_xlim(5, 50)
    y_max = max(max(df_ref['r']), max(df_alt['r']))
    # ax.set_ylim(0, y_max * 1.05)
    ax.grid(True)

    ax.plot([0, y_max], [0, y_max], c='grey', label='y=x', linestyle='dashed')

    x = [row['r'] for index, row in df_alt.iterrows() if index in df_ref.index]
    y = [row['r'] for index, row in df_ref.iterrows() if index in df_alt.index]

    ax.scatter(x=x, y=y, color='C1')

    ax.set_xlabel('Fitted r value for fixed Alt')
    ax.set_ylabel('Fitted r value for fixed Ref')

    plt.title('BAD={}'.format(BAD_dict[BAD]))

    slope, intercept, r_value, p_value, std_err = st.linregress(x, y)
    print(slope, intercept)
    ax.plot(x, np.array(x) * slope + intercept, color='#DC3220', label='y={:.2f} x {} {:.2f}'
            .format(slope, '+' if np.sign(intercept) > 0 else '-', abs(intercept)))

    plt.legend()

    plt.savefig(os.path.expanduser('~/AC_10/Figure_AS_10_r_scatter_ref_alt_{:.2f}.svg'.format(BAD)))
    plt.close(fig)
