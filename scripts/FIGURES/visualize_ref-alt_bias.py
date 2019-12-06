import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

stats = pd.read_table('bias_statistics.tsv')
#stats = stats[stats['allele_reads'] <= 1000]

fig, ax = plt.subplots(figsize=(10, 8))

sns.scatterplot(x=stats['allele_reads'], y=np.log(stats['ref']), color='blue', ax=ax)
sns.scatterplot(x=stats['allele_reads'], y=np.log(stats['alt']), color='green', ax=ax)
#sns.scatterplot(x=stats['allele_reads'], y=np.log10(stats['ref']-stats['alt']), ax=ax)

plt.grid(True)
plt.title('ref-alt bias for CELLLINE')
plt.ylabel('log10 count')
plt.savefig('bias_CELLLINE.png')
plt.show()