import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

ref_stats = pd.read_table('ref_statistics.tsv')
alt_stats = pd.read_table('alt_statistics.tsv')

fig, ax = plt.subplots(figsize=(10, 8))

sns.scatterplot(x=ref_stats['ref'], y=np.log10(ref_stats['count']), color='blue', ax=ax)
sns.scatterplot(x=alt_stats['alt'], y=np.log10(alt_stats['count']), color='green', ax=ax)

plt.grid(True)
plt.show()
