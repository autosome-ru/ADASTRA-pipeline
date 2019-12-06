import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

stats = pd.read_table('bias_statistics.tsv')

fig, ax = plt.subplots(figsize=(10, 8))

sns.scatterplot(x=stats['allele_reads'], y=np.log10(stats['ref']), color='blue', ax=ax)
sns.scatterplot(x=stats['allele_reads'], y=np.log10(stats['alt']), color='green', ax=ax)

plt.grid(True)
plt.show()
