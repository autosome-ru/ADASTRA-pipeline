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

    agr = True
    mode = 'by'

    if agr:
        stats = pd.read_table(os.path.expanduser('~/fdr_pvalue_bias_statistics_{}.tsv'.format(mode)))

        print(stats)

        fig, ax = plt.subplots(figsize=(10, 8))

        x = np.linspace(0, 1, 50)

        palts = stats.groupby('alt_p', as_index=False)['counts'].sum()
        # sns.barplot(x=x[:-1],
        #             y=[palts[(x[i] <= palts['alt_p']) & (palts['alt_p'] <= x[i + 1])]['counts'].sum()
        #                for i in range(len(x) - 1)], label='alt', ax=ax, color='C1', alpha=0.5)

        prefs = stats.groupby('ref_p', as_index=False)['counts'].sum()
        # sns.barplot(x=x[:-1],
        #             y=[prefs[(x[i] <= prefs['ref_p']) & (prefs['ref_p'] <= x[i + 1])]['counts'].sum()
        #                for i in range(len(x) - 1)], label='ref', ax=ax, color='C0', alpha=0.5)

        sns.barplot(x=x[:-1],
                    y=[prefs[(x[i] <= prefs['ref_p']) & (prefs['ref_p'] < x[i + 1])]['counts'].sum() +
                       palts[(x[i] <= palts['alt_p']) & (palts['alt_p'] < x[i + 1])]['counts'].sum()
                       for i in range(len(x) - 1)], label='sum', ax=ax, color='grey')

        ax.xaxis.set_major_locator(ticker.FixedLocator(np.arange(0, len(x), 5)))
        ax.xaxis.set_major_formatter(ticker.FixedFormatter(['{:.1f}'.format(f) for f in x[::5]]))
        ax.tick_params(axis="x", rotation=90)

        plt.grid(True)
        plt.legend()
        plt.title('p-value distribution agregated mode={}'.format(mode))
        plt.xlabel('p-value')
        plt.ylabel('snp count')

        plt.savefig(os.path.expanduser('~/fixed_alt/agr_pvalue_dist_mode={}.png'.format(mode)))

    else:
        for BAD in states:
            stats = pd.read_table(os.path.expanduser('~/pvalue_bias_statistics_BAD={:.1f}.tsv'.format(BAD)))

            fig, ax = plt.subplots(figsize=(10, 8))

            x = np.linspace(0, 1, 50)

            palts = stats.groupby('alt_p', as_index=False)['counts'].sum()
            # sns.barplot(x=x[:-1],
            #             y=[palts[(x[i] <= palts['alt_p']) & (palts['alt_p'] <= x[i + 1])]['counts'].sum()
            #                for i in range(len(x) - 1)], label='alt', ax=ax, color='C1', alpha=0.5)

            prefs = stats.groupby('ref_p', as_index=False)['counts'].sum()
            # sns.barplot(x=x[:-1],
            #             y=[prefs[(x[i] <= prefs['ref_p']) & (prefs['ref_p'] <= x[i + 1])]['counts'].sum()
            #                for i in range(len(x) - 1)], label='ref', ax=ax, color='C0', alpha=0.5)

            sns.barplot(x=x[:-1],
                        y=[prefs[(x[i] < prefs['ref_p']) & (prefs['ref_p'] <= x[i + 1])]['counts'].sum() +
                           palts[(x[i] < palts['alt_p']) & (palts['alt_p'] <= x[i + 1])]['counts'].sum()
                           for i in range(len(x) - 1)], label='sum', ax=ax, color='grey')

            ax.xaxis.set_major_locator(ticker.FixedLocator(np.arange(0, len(x), 5)))
            ax.xaxis.set_major_formatter(ticker.FixedFormatter(['{:.1f}'.format(f) for f in x[::5]]))
            ax.tick_params(axis="x", rotation=90)

            plt.grid(True)
            plt.legend()
            plt.title('p-value distribution for BAD={:.1f}'.format(BAD))
            plt.xlabel('p-value')
            plt.ylabel('snp count')

            plt.savefig(os.path.expanduser('~/fixed_alt/pvalue_dist_BAD={:.1f}.png'.format(BAD)))
