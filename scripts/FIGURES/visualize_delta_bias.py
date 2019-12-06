import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

#n = int(sys.argv[1])
n = 75
stats = pd.read_table('~/cover_bias_statistics.tsv')
stats = stats[stats['cover'] == n]

stats['ref'] = (stats['cover'] + stats['delta']) / 2
sns.barplot(x='ref', y='counts', data=stats)

plt.grid(True)
plt.title('ref-alt bias for cell line n={}'.format(n))
plt.ylabel('count')
plt.xlabel('ref_read_counts')
plt.xticks([])
#plt.savefig('ref-alt_bias_n={}.png'.format(n))
plt.show()
