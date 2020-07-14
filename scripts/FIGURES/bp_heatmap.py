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
idx = 0

for column in df1.columns:
    print(df1[column])
    df1[column] = df1[column].apply(lambda x: int(x.split(',')[idx]))

df2 = pd.read_table(os.path.expanduser('~/statistics_for_TFs_{}.tsv'.format('alt')))
df2.index = df2.iloc[:, 0]
df2 = df2.drop(df2.columns[0], 1)
for column in df2.columns:

    df2[column] = df2[column].apply(lambda x: int(x.split(',')[idx]))

df = df1 + df2
print(df)
# df = df.astype(float)
sns.heatmap(df, annot=True, fmt='d')
mydict = {1: 'all SNPs', 0: 'ASB SNPs'}
plt.ylabel('N datasets')
plt.xlabel('Max cover tr')

plt.grid(True)
plt.title('All TFs, database size {}'.format(mydict[idx]))
plt.savefig(os.path.expanduser('~/All_TFs_heat.png'))
plt.show()
