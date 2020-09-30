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
plt.rcParams["legend.framealpha"] = 0.5
plt.rcParams['axes.xmargin'] = 0
plt.rcParams['axes.ymargin'] = 0
plt.rcParams["legend.framealpha"] = 1


# PARAMS
lw = 0.05
a = 0.3
a1 = 0.8
states = []
mults = []
k562 = []
mcf7 = []
a549 = []
other = []
for sig in ['all', 'all_5', '123456', 'all_but_1.33_2.5']:
    for mult in [2, 3, 4, 5, 6, 7, 8, 9, 10]:
        print(sig, mult)
        model = "cor_by_snp_CAIC@{}@{:.1f}".format(sig, mult)

        df = pd.read_table(os.path.expanduser("~/Desktop/soos_BAD/cor_stats_test.tsv"))
        # df = df[df['total_snps'] >= 10000]
        df['color'] = df.apply(get_color, axis=1)
        df['hue'] = df.apply(get_hue, axis=1)
        df = df.sort_values('color', axis=0)
        df = df.fillna(value=0.01)
        df = df.dropna(subset=[model])
        df = df.sort_values(model, axis=0, ascending=False)
        df["delta_tau"] = df[model] - df["cor_by_snp_probe_CGH"]

        df_k562 = df.apply(lambda x: get_cl(x, "K562", model), axis=1)
        df_mcf7 = df.apply(lambda x: get_cl(x, "MCF7", model), axis=1)
        df_a549 = df.apply(lambda x: get_cl(x, "A549", model), axis=1)
        df_other = df.apply(lambda x: get_cl(x, "other", model), axis=1)

        # Draw 3 colors for cell lines vs cosmic
        fig, ax = plt.subplots()
        fig.tight_layout(pad=2)
        # ax.margins(x=0, y=0)

        x_k562, y_k562 = get_rect(df_k562, model)
        x_mcf7, y_mcf7 = get_rect(df_mcf7, model)
        x_other, y_other = get_rect(df_other, model)
        x_a549, y_a549 = get_rect(df_a549, model)

        ax.stackplot(x_k562, y_k562, alpha=a1,
                     linewidth=lw, color='C1', labels=["K562"])
        ax.stackplot(x_mcf7, y_mcf7, alpha=a1,
                     linewidth=lw, color='C2', labels=["MCF7"])
        ax.stackplot(x_other, y_other, alpha=a1,
                     linewidth=lw, color='C0', labels=["Other"])
        ax.stackplot(x_a549, y_a549, alpha=a1,
                     linewidth=lw, color='C3', labels=["A549"])

        ax.grid(True)
        ax.legend()
        ax.set_ylabel("Kendall's τ (Segmentation, COSMIC)")
        ax.set_xlabel("Dataset groups sorted by τ")
        ax.set_ylim(-0.4, 1)

        plt.title('{}@{:.1f}'.format(sig, mult))

        plt.savefig(os.path.expanduser('~/Desktop/soos_BAD/stack_plot@{}@{:.1f}.png'.format(sig, mult)), dpi=300)
        plt.close(fig)

        # Draw 3 colors for cell lines vs cosmic separately
        fig, ax = plt.subplots()
        fig.tight_layout(pad=2)
        # ax.margins(x=0, y=0)

        ax.hist(df_k562[(df_k562[model] != 0) & (df_k562['total_snps'] > 1000)][model], alpha=a, linewidth=lw,  color='C1', label="K562", bins=50, range=(-0.4, 1))
        ax.hist(df_mcf7[(df_mcf7[model] != 0) & (df_mcf7['total_snps'] > 1000)][model], alpha=a, linewidth=lw, color='C2', label="MCF7", bins=50, range=(-0.4, 1))
        ax.hist(df_other[(df_other[model] != 0) & (df_other['total_snps'] > 1000)][model], alpha=a, linewidth=lw, color='C0', label="Other", bins=50, range=(-0.4, 1))
        ax.hist(df_a549[(df_a549[model] != 0) & (df_a549['total_snps'] > 1000)][model], alpha=a, linewidth=lw, color='C3', label="A549", bins=50, range=(-0.4, 1))

        print(model)
        print('K562: {:.3f}'.format(np.mean(df_k562[(df_k562[model] != 0) & (df_k562[model] != 0.01)][model])))
        print('MCF7: {:.3f}'.format(np.mean(df_mcf7[(df_mcf7[model] != 0) & (df_mcf7[model] != 0.01)][model])))
        print('A549: {:.3f}'.format(np.mean(df_a549[(df_a549[model] != 0) & (df_a549[model] != 0.01)][model])))
        print('Other: {:.3f}'.format(np.mean(df_other[(df_other[model] != 0) & (df_other[model] != 0.01)][model])))

        states.append(sig)
        mults.append(mult)
        k562.append(np.quantile(df_k562[(df_k562[model] != 0) & (df_k562[model] != 0.01)][model], 0.75))
        mcf7.append(np.quantile(df_mcf7[(df_mcf7[model] != 0) & (df_mcf7[model] != 0.01)][model], 0.75))
        a549.append(np.quantile(df_a549[(df_a549[model] != 0) & (df_a549[model] != 0.01)][model], 0.75))
        other.append(np.quantile(df_other[(df_other[model] != 0) & (df_other[model] != 0.01)][model], 0.75))

        ax.grid(True)
        ax.legend()
        ax.set_ylabel("Dataset groups count")
        ax.set_xlabel("Kendall's τ (Segmentation, COSMIC)")

        plt.title('{}@{:.1f}'.format(sig, mult))

        plt.savefig(os.path.expanduser('~/Desktop/soos_BAD/hist_plot@{}@{:.1f}.png'.format(sig, mult)), dpi=300)
        plt.close(fig)

        # Draw scatter vs COSMIC
        fig, ax = plt.subplots()
        fig.tight_layout(pad=2)

        df_k562 = df[df['color'] == 'C1']
        df_mcf7 = df[df['color'] == 'C2']
        df_a549 = df[df['color'] == 'C3']
        df_other = df[df['color'] == 'C0']

        field = 'total_snps'

        sc1 = plt.scatter(y=df_k562[model], x=df_k562[field], zorder=1,
                        linewidth=0, alpha=0.7, color='C1', label='K562')
        sc2 = plt.scatter(y=df_mcf7[model], x=df_mcf7[field], zorder=1,
                        linewidth=0, alpha=0.7, color='C2', label='MCF7')
        sc3 = plt.scatter(y=df_other[model], x=df_other[field], zorder=1,
                        linewidth=0, alpha=0.7, color='C0', label='Other')
        sc4 = plt.scatter(y=df_a549[model], x=df_a549[field], zorder=1,
                        linewidth=0, alpha=0.7, color='C3', label='A 549')
        sns.lineplot(x=[-1, 1], y=[-1, 1], color='#505050')
        ax.axvline(x=0, color='#505050', linestyle='--')
        ax.axhline(y=0, color='#505050', linestyle='--')
        ax.legend(loc='lower right', handletextpad=0.3, handlelength=1)
        # ax.set_xticks([-0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])
        ax.set_yticks([-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])
        ax.set_xscale('log')
        ax.set_xlim(100, 2000000)
        ax.set_ylim(-0.5, 1)
        ax.grid(True)

        ax.set_ylabel("Kendall's τ (Segmentation, COSMIC)")
        ax.set_xlabel("Number of SNPs in a group of datasets")

        plt.title('{}@{:.1f}'.format(sig, mult))

        plt.savefig(os.path.expanduser('~/Desktop/soos_BAD/cor_plot@{}@{:.1f}.png'.format(sig, mult)), dpi=300)
        plt.close(fig)

df = pd.DataFrame({
    'States': states,
    'Multiplier': mults,
    'K562_Q3': k562,
    'MCF7_Q3': mcf7,
    'A549_Q3': a549,
    'Other': other,
})

df.to_csv(os.path.expanduser('~/Desktop/soos_BAD/stats.tsv'), sep='\t', index=False)
