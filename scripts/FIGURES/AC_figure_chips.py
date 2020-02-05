import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, ticker
import seaborn as sns


def get_color(row):
    print(row)
    if row['#cell_line'] == 'K562__myelogenous_leukemia_':
        return 'K562__myelogenous_leukemia_'
    elif row['#cell_line'] == "MCF7__Invasive_ductal_breast_carcinoma_":
        return "MCF7__Invasive_ductal_breast_carcinoma_"
    else:
        return "others"


df = pd.read_table(os.path.expanduser("~/cor_stats_test.tsv"))
df = df[df['total_snps'] >= 30000]
df['color'] = df.apply(get_color, axis=1)
df = df.sort_values('color', axis=0)
sns.scatterplot(x="cor_by_snp_CAIC", y="cor_by_snp_probe_CGH", hue='color', data=df,  linewidth=0, alpha=0.5)
sns.lineplot(x=[-1, 1], y=[-1, 1], color='#505050')
plt.grid(True)
plt.xlim(-1, 1)
plt.ylim(-1, 1)
plt.savefig(os.path.expanduser("~/AC_chips.png"), dpi=300)
plt.close()
