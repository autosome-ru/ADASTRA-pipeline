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
        return 'C1'
    elif row['#cell_line'] == "MCF7__Invasive_ductal_breast_carcinoma_":
        return "C2"
    else:
        return "C0"


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


def get_rect(dataframe: pd.DataFrame, col_name):
    tri_list = dataframe[col_name].to_list()
    x = [0]
    list = [tri_list[0]]
    for i in range(1, len(tri_list)):

        if tri_list[i] != 0 and tri_list[i-1] == 0:
            list.append(0)
            x.append(i - 0.5)
            list.append(tri_list[i])
            x.append(i - 0.5)

        elif tri_list[i] == 0 and tri_list[i-1] != 0:
            list.append(tri_list[i - 1])
            x.append(i - 0.5)
            list.append(0)
            x.append(i - 0.5)
        else:
            list.append(tri_list[i])
            x.append(i)
    return x, list


sns.set(font_scale=1.4, style="ticks", font="lato",
        # palette=('#56B4E9', '#009E73', '#F0E442'))
        # palette=('#7570b3', '#d95f02', '#1b9e77'))
        palette=('#56B4E9', '#E69F00', '#009E73'))
        # palette=('#1f77b4', '#2ca02c', '#ff7f0e'))
sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})
plt.rcParams['font.weight'] = "medium"
plt.rcParams['axes.labelweight'] = 'medium'
plt.rcParams['figure.titleweight'] = 'medium'
plt.rcParams['axes.titleweight'] = 'medium'
plt.rcParams['figure.figsize'] = 6, 5
plt.rcParams["legend.framealpha"] = 1
plt.rcParams['axes.xmargin'] = 0
plt.rcParams['axes.ymargin'] = 0
plt.rcParams["legend.framealpha"] = 1


# PARAMS
lw = 0.05
a = 0.8


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
fig.tight_layout(pad=2)
# ax.margins(x=0, y=0)

x_k562, y_k562 = get_rect(df_k562, "cor_by_snp_CAIC")
x_mcf7, y_mcf7 = get_rect(df_mcf7, "cor_by_snp_CAIC")
x_other, y_other = get_rect(df_other, "cor_by_snp_CAIC")

ax.stackplot(x_k562, y_k562, alpha=a,
             linewidth=lw, color='C1', labels=["K562"])
ax.stackplot(x_mcf7, y_mcf7, alpha=a,
             linewidth=lw, color='C2', labels=["MCF7"])
ax.stackplot(x_other, y_other, alpha=a,
             linewidth=lw, color='C0', labels=["Other"])

ax.grid(True)
ax.legend()
ax.set_ylabel("Kendall's τ (Segmentation, COSMIC)")
ax.set_xlabel("Dataset groups sorted by τ")

plt.savefig(os.path.expanduser("~/AC_9/Figure_AS_9_cor_cosmic.png"), dpi=300)
plt.savefig(os.path.expanduser("~/AC_9/Figure_AS_9_cor_cosmic.svg"), dpi=300)
plt.close(fig)

# Draw scatter vs COSMIC
fig, ax = plt.subplots()
fig.tight_layout(pad=2)

sns.scatterplot(x="cor_by_snp_CAIC", y="cor_by_snp_probe_CGH", zorder=10,
                data=df[df['color'] == 'C1'], linewidth=0, alpha=0.7, color='C1', label='K562')
sns.scatterplot(x="cor_by_snp_CAIC", y="cor_by_snp_probe_CGH", zorder=10,
                data=df[df['color'] == 'C2'], linewidth=0, alpha=0.7, color='C2', label='MCF7')
sns.scatterplot(x="cor_by_snp_CAIC", y="cor_by_snp_probe_CGH", zorder=10,
                data=df[df['color'] == 'C0'], linewidth=0, alpha=0.7, color='C0', label='Other')
sns.lineplot(x=[-1, 1], y=[-1, 1], color='#505050')
ax.axvline(x=0, color='#505050', linestyle='--')
ax.axhline(y=0, color='#505050', linestyle='--')
ax.legend(loc='lower left')
ax.set_xticks([-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])
ax.set_yticks([-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])
ax.set_xlim(-0.5, 1)
ax.set_ylim(-0.5, 1)
ax.grid(True)

ax.set_ylabel("Kendall's τ (aCGH Varma et al., COSMIC)")
ax.set_xlabel("Kendall's τ (Segmentation, COSMIC)")

plt.savefig(os.path.expanduser("~/AC_9/Figure_AS_9_scatter.png"), dpi=300)
plt.savefig(os.path.expanduser("~/AC_9/Figure_AS_9_scatter.svg"), dpi=300)
plt.close(fig)

# delta tau 3 colors
fig, ax = plt.subplots()
fig.tight_layout(pad=2)

df = df.sort_values("delta_tau", axis=0, ascending=False)
df = df.dropna(subset=["cor_by_snp_probe_CGH"])

df_k562 = df.apply(lambda x: get_cl(x, "K562"), axis=1)
df_mcf7 = df.apply(lambda x: get_cl(x, "MCF7"), axis=1)
df_other = df.apply(lambda x: get_cl(x, "other"), axis=1)
x_k562, y_k562 = get_rect(df_k562, "delta_tau")
x_mcf7, y_mcf7 = get_rect(df_mcf7, "delta_tau")
x_other, y_other = get_rect(df_other, "delta_tau")

ax.stackplot(x_k562, y_k562, alpha=a,
             linewidth=lw, color='C1', labels=["K562"])
ax.stackplot(x_mcf7, y_mcf7, alpha=a,
             linewidth=lw, color='C2', labels=["MCF7"])
ax.stackplot(x_other, y_other, alpha=a,
             linewidth=lw, color='C0', labels=["Other"])

ax.set_ylabel("Δτ (Segmentation, aCGH Varma et al.)")
ax.set_xlabel("Dataset groups sorted by Δτ")

ax.grid(True)
ax.legend()
plt.savefig(os.path.expanduser("~/AC_9/Figure_AS_9_delta_tau_chips.png"), dpi=300)
plt.savefig(os.path.expanduser("~/AC_9/Figure_AS_9_delta_tau_chips.svg"), dpi=300)
plt.close(fig)
