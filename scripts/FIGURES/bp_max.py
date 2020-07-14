import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_style("white")
sns.set_context({"figure.figsize": (10, 8)})

df1 = pd.read_table('statistics_for_TFs_{}_max_cov.tsv'.format('ref'))
df1 = df1.drop('#name',1)
df1 = df1.drop('name',1)

df2 = pd.read_table('statistics_for_TFs_{}_max_cov.tsv'.format('alt'))
df2 = df2.drop('#name',1)
df2 = df2.drop('name',1)

sum_n1 = pd.DataFrame(df1.sum(axis=0).transpose())
sum_n1["tr"] = sum_n1.index.astype(int)
sum_n1.columns = ["counts", "tr"]
print(sum_n1.head())

sum_n2 = pd.DataFrame(df2.sum(axis=0).transpose())
sum_n2["tr"] = sum_n2.index.astype(int)
sum_n2.columns = ["counts", "tr"]
print(sum_n2.head())


sns.barplot(x=sum_n1.tr, y=sum_n2.counts + sum_n1.counts, color = "C1")
bottom_plot = sns.barplot(x=sum_n1.tr, y=sum_n1.counts, color = "C0")


topbar = plt.Rectangle((0,0),1,1,fc="C1", edgecolor = 'none')
bottombar = plt.Rectangle((0,0),1,1,fc='C0',  edgecolor = 'none')
l = plt.legend([bottombar, topbar], ['Ref', 'Alt'], loc=1, ncol = 2, prop={'size':16})
l.draw_frame(False)

bottom_plot.set_ylabel("SNP FDR <= 0.05")
bottom_plot.set_xlabel("x: max_cover >= x")

plt.grid(True)
plt.title('All TFs, database size')
plt.savefig('All_TFs_max.png')
plt.show()
