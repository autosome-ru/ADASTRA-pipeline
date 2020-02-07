import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, ticker
import seaborn as sns


def get_hue(row):
    if row['#cell_line'] == 'K562__myelogenous_leukemia_':
        return 'K562__myelogenous_leukemia_'
    elif row['#cell_line'] == "MCF7__Invasive_ductal_breast_carcinoma_":
        return "MCF7__Invasive_ductal_breast_carcinoma_"
    else:
        return "others"


def get_color(row):
    if row['#cell_line'] == 'K562__myelogenous_leukemia_':
        return 'C0'
    elif row['#cell_line'] == "MCF7__Invasive_ductal_breast_carcinoma_":
        return "C1"
    else:
        return "C2"


def get_cl(row, name):
    if name == "K562":
        if row['#cell_line'] != 'K562__myelogenous_leukemia_':
            row["cor_by_snp_CAIC"] = 0
            row["delta_tau"] = 0
    if name == "MCF7":
        if row['#cell_line'] != "MCF7__Invasive_ductal_breast_carcinoma_":
            row["cor_by_snp_CAIC"] = 0
            row["delta_tau"] = 0
    elif name == "other":
        if row['#cell_line'] == "MCF7__Invasive_ductal_breast_carcinoma_" or \
                row['#cell_line'] == 'K562__myelogenous_leukemia_':
            row["cor_by_snp_CAIC"] = 0
            row["delta_tau"] = 0
    return row


df = pd.read_table(os.path.expanduser("~/cor_stats_test.tsv"))
# df = df[df['total_snps'] >= 10000]
df['color'] = df.apply(get_color, axis=1)
df['hue'] = df.apply(get_hue, axis=1)
df = df.sort_values('color', axis=0)
df = df.dropna(subset=["cor_by_snp_CAIC"])
df = df.sort_values("cor_by_snp_CAIC", axis=0, ascending=False)
df["delta_tau"] = df["cor_by_snp_CAIC"] - df["cor_by_snp_probe_CGH"]

df_k562 = df.apply(lambda x: get_cl(x, "K562"), axis=1)
df_mcf7 = df.apply(lambda x: get_cl(x, "MCF7"), axis=1)
df_other = df.apply(lambda x: get_cl(x, "other"), axis=1)

# Draw 3 colors for cell lines vs cosmic
fig, ax = plt.subplots()
fig.tight_layout(pad=1.5)
# ax.margins(x=0, y=0)

ax.stackplot(range(len(df.index)), df_k562["cor_by_snp_CAIC"], alpha=0.7,
             linewidth=0.7, color='C1', labels=["K562"])
ax.stackplot(range(len(df.index)), df_mcf7["cor_by_snp_CAIC"], alpha=0.7,
             linewidth=0.7, color='C2', labels=["MCF7"])
ax.stackplot(range(len(df.index)), df_other["cor_by_snp_CAIC"], alpha=0.7,
             linewidth=0.7, color='C0', labels=["other"])

ax.grid(True)
ax.legend()
plt.savefig(os.path.expanduser("~/AC_9/AS_Figure_9_cor_cosmic.png"), dpi=300)
plt.close(fig)

# Draw scatter vs COSMIC
fig, ax = plt.subplots()
fig.tight_layout(pad=1.5)

sns.scatterplot(x="cor_by_snp_CAIC", y="cor_by_snp_probe_CGH", hue='hue', data=df, linewidth=0, alpha=0.5)
sns.lineplot(x=[-1, 1], y=[-1, 1], color='#505050')
ax.axvline(x=0, ymax=0, color='#505050', linestyle='--')
ax.axhline(y=0, xmax=0, color='#505050', linestyle='--')
ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)
ax.grid(True)
plt.savefig(os.path.expanduser("~/AC_9/AS_Figure_9_scatter.png"), dpi=300)
plt.close(fig)


# delta tau 3 colors
fig, ax = plt.subplots()
fig.tight_layout(pad=1.5)

df = df.sort_values("delta_tau", axis=0, ascending=False)
df = df.dropna(subset=["cor_by_snp_probe_CGH"])

df_k562 = df.apply(lambda x: get_cl(x, "K562"), axis=1)
df_mcf7 = df.apply(lambda x: get_cl(x, "MCF7"), axis=1)
df_other = df.apply(lambda x: get_cl(x, "other"), axis=1)

ax.stackplot(range(len(df.index)), df_k562["delta_tau"], alpha=0.7,
              linewidth=0.7, color='C1', labels=["K562"])
ax.stackplot(range(len(df.index)), df_mcf7["delta_tau"], alpha=0.7,
              linewidth=0.7, color='C2', labels=["MCF7"])
ax.stackplot(range(len(df.index)), df_other["delta_tau"], alpha=0.7,
              linewidth=0.7, color='C0', labels=["other"])

ax.grid(True)
ax.legend()
plt.savefig(os.path.expanduser("~/AC_9/AS_Figure_9_delta_tau_chips.png"), dpi=300)
plt.close(fig)

