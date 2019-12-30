import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

for BAD in [1, 2, 3, 4, 5, 6, 4/3, 5/2, 3/2]:
    stats = pd.read_table(os.path.expanduser('~/bias_statistics_BAD={}.tsv'.format(BAD)))
    # stats = stats[stats['allele_reads'] <= 1000]
    fig, ax = plt.subplots(figsize=(10, 8))

    sns.scatterplot(x=stats['allele_reads'], y=np.log(stats['ref']), color='blue', ax=ax)
    sns.scatterplot(x=stats['allele_reads'], y=np.log(stats['alt']), color='green', ax=ax)
    #sns.scatterplot(x=stats['allele_reads'], y=np.log10(stats['ref']-stats['alt']), ax=ax)

    plt.grid(True)
    plt.title('ref-alt bias for CELLLINE')
    plt.ylabel('log10 count')
    plt.savefig('bias_CELLLINE.png')
    plt.show()
