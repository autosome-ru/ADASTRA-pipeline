import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

max_c = 50

t = pd.read_table(os.path.expanduser('~/diploid_snps_statistics.tsv'))
t = t[(t['ref'] <= max_c) & (t['alt'] <= max_c)]
t = t.pivot('ref', 'alt', 'count')
t.sort_index(ascending=False, inplace=True)
t.fillna(0, inplace=True)

for k in range(max_c + 1):
    for l in range(k + 1):
        t[k][l], t[l][k] = t[k][l] / (t[k][l] + t[l][k]), t[l][k] / (t[k][l] + t[l][k])

sns.heatmap(t)
plt.show()
