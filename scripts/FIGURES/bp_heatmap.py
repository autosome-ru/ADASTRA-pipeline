import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_style("white")
sns.set_context({"figure.figsize": (10, 8)})

df1 = pd.read_table(os.path.expanduser('~/statistics_for_TFs_{}.tsv'.format('ref')))
df1.index = df1.iloc[:, 0]
df1 = df1.drop(df1.columns[0], 1)

df2 = pd.read_table(os.path.expanduser('~/statistics_for_TFs_{}.tsv'.format('alt')))
df2.index = df2.iloc[:, 0]
df2 = df2.drop(df2.columns[0], 1)

df = df1 + df2
df = df.astype(int)
sns.heatmap(df, annot=True, fmt='d')

plt.ylabel('Total cover tr')
plt.xlabel('Max cover tr')

plt.grid(True)
plt.title('All TFs, database size')
plt.savefig(os.path.expanduser('~/All_TFs_heat.png'))
plt.show()
