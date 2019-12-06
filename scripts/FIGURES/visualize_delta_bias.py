import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

n = int(sys.argv[1])
stats = pd.read_table('cover_bias_statistics.tsv')
stats = stats[stats['cover'] == n]
stats['ref'] = (stats['counts'] + stats['delta']) / 2
sns.barplot(stats[['ref', 'counts']])

plt.grid(True)
plt.title('ref-alt bias for CELLLINE')
plt.ylabel('count')
plt.xlabel('delta(ref-alt)')
plt.savefig('bias_CELLLINE.png')
plt.show()
