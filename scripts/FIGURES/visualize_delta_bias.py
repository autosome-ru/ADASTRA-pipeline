import sys
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

#n = int(sys.argv[1])

for n in range(10, 101, 5):
    stats = pd.read_table('~/cover_bias_statistics.tsv')
    stats = stats[stats['cover'] == n]

    sns.barplot(x='ref', y='counts', data=stats)

    plt.grid(True)
    plt.title('ref-alt bias for H1_embryonic_stem_cells n={}'.format(n))
    plt.ylabel('count')
    plt.xlabel('ref_read_counts')
    plt.savefig('/home/serj/ref-alt_bias_n-{}.png'.format(n))
plt.show()
