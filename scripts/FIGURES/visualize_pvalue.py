import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

stats = pd.read_table(os.path.expanduser('~/pvalue_bias_statistics_BAD=1.0.tsv'))
print(stats)

fig, ax = plt.subplots(figsize=(10, 8))

prefs = stats[stats['ref_p'] != 1]
plt.plot(x=prefs['ref_p'], y=prefs['counts'], label='ref')

palts = stats[stats['alt_p'] != 1]
plt.plot(x=palts['alt_p'], y=palts['counts'], label='alt')

plt.grid(True)
plt.title('ref-alt p_value')
plt.xlabel('p_value')
plt.ylabel('count')
plt.savefig(os.path.expanduser('~/fixed_alt/p_dist.png'))
plt.show()
