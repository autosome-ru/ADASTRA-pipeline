import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_style("white")
sns.set_context({"figure.figsize": (10, 8)})

df1 = pd.read_table(os.path.expanduser('~/statistics_for_TFs_weight_ref.tsv'))
df1 = df1.drop("Unnamed: 0", 1)
print(df1)



sns.barplot(x=df[0], y=df[1], color="C1")
bottom_plot = sns.barplot(x=df[0], y=df[1], color = "C0")


topbar = plt.Rectangle((0,0),1,1,fc="C1", edgecolor = 'none')
bottombar = plt.Rectangle((0,0),1,1,fc='C0',  edgecolor = 'none')
l = plt.legend([bottombar, topbar], ['Ref', 'Alt'], loc=1, ncol = 2, prop={'size':16})
l.draw_frame(False)

bottom_plot.set_ylabel("SNP FDR <= 0.05")
bottom_plot.set_xlabel("x: min_weight >= x")

plt.grid(True)
plt.title('All TFs, database size')
plt.savefig('All_TFs_total.png')
plt.show()
