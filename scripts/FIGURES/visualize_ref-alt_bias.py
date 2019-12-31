import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


def quantileNormalize(df_input):
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df


def quantileNormalizeToAlt(df_input):
    df = df_input.copy()
    #compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df['alt'].tolist()
    #sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df


flag = 'after'

#for flag in ('before', 'after', 'ref', 'alt'):
for flag in ('ref', 'alt'):
    for BAD in [1, 2, 3, 4, 5, 6, 4/3, 5/2, 3/2]:
        stats = pd.read_table(os.path.expanduser('~/bias_statistics_BAD={:.1f}.tsv'.format(BAD)))
        stats = stats[stats['allele_reads'] >= 500]
        fig, ax = plt.subplots(figsize=(10, 8))

        df_tonorm = stats[['ref', 'alt']]
        new_stats = pd.DataFrame()
        new_stats[['ref', 'alt']] = quantileNormalizeToAlt(df_tonorm)

        if flag in ('ref', 'before'):
            sns.scatterplot(x=stats['allele_reads'], y=stats['ref'], color='C0', ax=ax, label='ref', s=50)
        if flag in ('ref', 'after'):
            sns.scatterplot(x=stats['allele_reads'], y=new_stats['ref'], color='C3', ax=ax, label='ref norm', s=40)
        if flag in ('alt', 'before'):
            sns.scatterplot(x=stats['allele_reads'], y=stats['alt'], color='C1', ax=ax, label='alt', s=40)
        if flag in ('alt', 'after'):
            sns.scatterplot(x=stats['allele_reads'], y=new_stats['alt'], color='C4', ax=ax, label='alt norm', s=40)
        if flag == 'ratio':
            sns.scatterplot(x=stats['allele_reads'], y=np.log10(stats['ref'] + 1)-np.log10(1 + stats['alt']), ax=ax, color='C5', label='log10 (ref+1)/(alt+1)', s=40)
            plt.axhline(y=0, color='black')

        plt.grid(True)
        if flag == 'ratio':
            plt.title('ref-alt ratio for BAD={:.1f}'.format(BAD))
            plt.ylabel('log10 (ref+1)/(alt+1)')
        else:
            plt.title('ref-alt bias for BAD={:.1f}'.format(BAD))
            plt.ylabel('log10 count')
        plt.legend()
        plt.savefig(os.path.expanduser('~/qnorm/bias_q_{}_BAD={:.1f}.png'.format(flag, BAD)))
        plt.close(fig)
