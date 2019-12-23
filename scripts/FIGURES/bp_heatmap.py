import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_style("white")
sns.set_context({"figure.figsize": (10, 8)})

df1 = pd.read_table('statistics_for_TFs_{}.tsv'.format('ref'))
df1.index = df1.columns

df2 = pd.read_table('statistics_for_TFs_{}.tsv'.format('alt'))
df2.index = df2.columns

df = df1 + df2
df = df.astype(int)
sns.heatmap(df, annot=True, fmt='d')

plt.ylabel('Total cover tr')
plt.xlabel('Max cover tr')

plt.grid(True)
plt.title('All TFs, database size')
plt.savefig('All_TFs_heat.png')
plt.show()
