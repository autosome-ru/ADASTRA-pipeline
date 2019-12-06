import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

stats = pd.read_table('bias_statistics.tsv')

fig, ax = plt.subplots(figsize=(10, 8))

#sns.scatterplot(x=stats['allele_reads'], y=np.log(stats['ref']), color='blue', ax=ax)
#sns.scatterplot(x=stats['allele_reads'], y=np.log(stats['alt']), color='green', ax=ax)
sns.scatterplot(x=stats['allele_reads'], y=(stats['ref']-stats['alt'])/np.sqrt(stats['alt']), ax=ax)

plt.grid(True)
plt.show()
