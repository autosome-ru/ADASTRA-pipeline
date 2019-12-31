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


for BAD in [1, 2, 3, 4, 5, 6, 4/3, 5/2, 3/2]:
    stats = pd.read_table(os.path.expanduser('~/bias_statistics_BAD={:.1f}.tsv'.format(BAD)))
    #stats = stats[stats['allele_reads'] <= 500]
    fig, ax = plt.subplots(figsize=(10, 8))

    print(stats[stats['ref'] == 0].size)

    df_tonorm = stats[['ref', 'alt']]
    stats[['ref', 'alt']] = quantileNormalize(df_tonorm)

    print(stats[stats['alt'] == 0.5].size)

    sns.scatterplot(x=stats['allele_reads'], y=np.log10(stats['ref']), color='blue', ax=ax, label='ref')
    sns.scatterplot(x=stats['allele_reads'], y=np.log10(stats['alt']), color='green', ax=ax, label='alt')
    #sns.scatterplot(x=stats['allele_reads'], y=np.log10(stats['ref']-stats['alt']), ax=ax)

    plt.grid(True)
    plt.title('ref-alt bias for BAD={:.1f}'.format(BAD))
    plt.ylabel('log10 count')
    plt.legend()
    plt.savefig(os.path.expanduser('~/qnorm/bias_q_BAD={:.1f}.png'.format(BAD)))
