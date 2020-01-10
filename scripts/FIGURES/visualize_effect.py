import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, ticker
import seaborn as sns

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import parameters_path
from scripts.HELPERS.helpers import states

if __name__ == '__main__':

    mode = 'maxdepth'

    stats1 = pd.read_table(os.path.expanduser('~/fdr_effect_size_gr_005_{}.tsv'.format(mode)))
    stats2 = pd.read_table(os.path.expanduser('~/fdr_effect_size_le_005_{}.tsv'.format(mode)))

    #stats1 = stats1[stats1['metric'] != 0]
    #stats2 = stats2[stats2['metric'] != 0]


    stats1['dens'] = stats1['counts']/stats1['counts'].sum()
    stats2['dens'] = stats2['counts'] / stats2['counts'].sum()

    print(stats1)
    #stats = stats1[stats1['metric'] < 500]
    fig, ax = plt.subplots(figsize=(10, 8))
    #stats1['counts'] = stats1['counts'].apply(lambda x: np.math.log10(x))
    # sns.scatterplot(x='metric', y='dens', data=stats1, ax=ax, label='NOASB')
    # sns.scatterplot(x='metric', y='dens', data=stats2, ax=ax, label='ASB')
    #
    x = np.linspace(-6, 6, 50)

    sns.barplot(x=x[:-1],
                y=[stats1[(x[i] <= stats1['metric']) & (stats1['metric'] <= x[i + 1])]['dens'].sum()
                   for i in range(len(x) - 1)], label='NOASB', ax=ax, color='C0', alpha=0.5)

    sns.barplot(x=x[:-1],
                y=[stats2[(x[i] <= stats2['metric']) & (stats2['metric'] <= x[i + 1])]['dens'].sum()
                   for i in range(len(x) - 1)], label='ASB', ax=ax, color='C1', alpha=0.5)

    # sns.barplot(x=x[:-1],
    #             y=[prefs[(x[i] <= prefs['ref_p']) & (prefs['ref_p'] < x[i + 1])]['counts'].sum() +
    #                palts[(x[i] <= palts['alt_p']) & (palts['alt_p'] < x[i + 1])]['counts'].sum()
    #                for i in range(len(x) - 1)], label='sum', ax=ax, color='grey')
    #
    ax.xaxis.set_major_locator(ticker.FixedLocator(np.arange(0, len(x), 5)))
    ax.xaxis.set_major_formatter(ticker.FixedFormatter(['{:.1f}'.format(f) for f in x[::5]]))
    ax.tick_params(axis="x", rotation=90)

    plt.grid(True)
    plt.legend()
    plt.title('{} effect size dist'.format(mode))
    plt.xlabel('effect size')
    plt.ylabel('snp count density')
    plt.savefig(os.path.expanduser('~/fixed_alt/effect_dist_{}.png'.format(mode)))
