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
    elif row['#cell_line'] == 'A549__lung_carcinoma_':
        return 'A549__lung_carcinoma_'
    else:
        return "others"


def get_color(row):
    if row['#cell_line'] == 'K562__myelogenous_leukemia_':
        return 'C1'
    elif row['#cell_line'] == "MCF7__Invasive_ductal_breast_carcinoma_":
        return "C2"
    elif row['#cell_line'] == 'A549__lung_carcinoma_':
        return 'C3'
    else:
        return "C0"


def get_cl(row, name, model):
    if name == "K562":
        if row['#cell_line'] != 'K562__myelogenous_leukemia_':
            row[model] = 0
            row["delta_tau"] = 0
    if name == "MCF7":
        if row['#cell_line'] != "MCF7__Invasive_ductal_breast_carcinoma_":
            row[model] = 0
            row["delta_tau"] = 0
    if name == 'A549':
        if row['#cell_line'] != 'A549__lung_carcinoma_':
            row[model] = 0
            row["delta_tau"] = 0
    elif name == "other":
        if row['#cell_line'] == "MCF7__Invasive_ductal_breast_carcinoma_" or \
                row['#cell_line'] == 'K562__myelogenous_leukemia_' or row['#cell_line'] == 'A549__lung_carcinoma_':
            row[model] = 0
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
        palette=('#56B4E9', '#E69F00', '#009E73', '#D55E00'))
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

# model = 'cor_by_snp_CAIC@{}@{:.1f}'.format('all_5', 4)
model = 'cor_by_snp_CAIC'

df = pd.read_table(os.path.expanduser("~/DataForFigures/cor_stats_test_full.tsv"))
# df = df[df['total_snps'] >= 10000]
df['color'] = df.apply(get_color, axis=1)
df['hue'] = df.apply(get_hue, axis=1)
df = df.sort_values('color', axis=0)
df = df.fillna(value=0.01)
df = df.dropna(subset=[model])
df = df.sort_values(model, axis=0, ascending=False)
df["delta_tau"] = df[model] - df["cor_by_snp_probe_CGH"]

# Draw scatter vs COSMIC
fig, ax = plt.subplots()
fig.tight_layout(pad=2)

field = 'total_snps'

sc = plt.scatter(y=df[model], x=df[field], zorder=1, linewidth=0, alpha=0.7, color=df['color'])

# df_k562 = df[df['color'] == 'C1']
# df_mcf7 = df[df['color'] == 'C2']
# df_a549 = df[df['color'] == 'C3']
# df_other = df[df['color'] == 'C0']
#
# field = 'total_snps'
#
# sc1 = plt.scatter(y=df_k562[model], x=df_k562[field], zorder=1,
#                   linewidth=0, alpha=0.7, color='C1', label='K562')
# sc2 = plt.scatter(y=df_mcf7[model], x=df_mcf7[field], zorder=1,
#                   linewidth=0, alpha=0.7, color='C2', label='MCF7')
# sc3 = plt.scatter(y=df_other[model], x=df_other[field], zorder=1,
#                   linewidth=0, alpha=0.7, color='C0', label='Other')
# sc4 = plt.scatter(y=df_a549[model], x=df_a549[field], zorder=1,
#                   linewidth=0, alpha=0.7, color='C3', label='A 549')



sns.lineplot(x=[-1, 1], y=[-1, 1], color='#505050')
ax.axvline(x=0, color='#505050', linestyle='--')
ax.axhline(y=0, color='#505050', linestyle='--')
ax.legend(loc='lower right', handletextpad=0.3, handlelength=1)
# ax.set_xticks([-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])
ax.set_yticks([-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])
ax.set_xscale('log')
ax.set_xlim(100, 4000000)
ax.set_ylim(-0.5, 1)
ax.grid(True)

ax.set_ylabel("Kendall's τ (Segmentation, COSMIC)")
ax.set_xlabel("Number of SNPs in a group of datasets")

annot = ax.annotate("", xy=(0, 0), xytext=(12, 12), textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->"), fontsize=10)
annot.set_visible(False)


def update_annot(ind, names, sc, col):
    pos = sc.get_offsets()[ind["ind"][0]]
    annot.xy = pos
    text = "{}".format("\n\n".join([names[n] for n in ind["ind"]]))
    annot.set_text(text)
    annot.get_bbox_patch().set_facecolor(col)
    annot.get_bbox_patch().set_alpha(0.4)


def hover_sc(names, sc, col):
    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = sc.contains(event)
            if cont:
                update_annot(ind, names, sc, col)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()
    return hover


def extract_biosamples(row):
    name = row['cells']
    if name.find('ENC') != -1:
        return name[name.find('ENC'):-1]
    if name.find('GSE') != -1:
        return name[name.find('GSE'):]


fig.canvas.mpl_connect("motion_notify_event", hover_sc(['\n'.join(map(str, [x['#cell_line'], x['#_of_merged_datasets'], extract_biosamples(x)])) for i, x in df.iterrows()], sc, 'grey'))

plt.show()
