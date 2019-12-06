import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

#n = int(sys.argv[1])

stats = pd.read_table('~/cover_bias_statistics.tsv')

for n in range(10, 101, 5):
    statsplot = pd.DataFrame(stats.loc[stats['cover'] == n])
    statsplot['ref_counts'] = statsplot['ref_counts'].astype(int)
    statsplot['counts'] = statsplot['counts'].astype(int)
    statsplot['color'] = np.abs(statsplot['ref_counts'] - n/2)
    print('made data for n={}'.format(n))
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.barplot(x='ref_counts', y='counts', data=statsplot, hue='color',
                ax=ax, dodge=False,
                palette='husl')

    plt.title('ref-alt bias for BAD=1 n={}'.format(n))
    ax.legend().remove()
    plt.ylabel('count')
    plt.xlabel('ref_read_counts')
    plt.savefig(os.path.expanduser('~/ref-alt_bias_BAD1_n-{}.png'.format(n)))
plt.show()
