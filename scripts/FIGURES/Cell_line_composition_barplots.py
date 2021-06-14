import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, ticker
import seaborn as sns


def get_state_label(string):
    f = float(string)
    if round(f) == f:
        return str(int(f))
    if round(2*f) == 2*f:
        return str(round(2*f))+'/2'
    if round(3*f + 0.01) == round(3*f + 0.01, 5) or round(3*(f-0.01) + 0.02) == round(3*(f-0.01) + 0.02, 5):
        return str(round(3*f))+'/3'
    if round(4*f) == 4*f:
        return str(round(4*f))+'/4'
    if round(5*f) == 5*f:
        return str(round(5*f))+'/5'
    return string


def get_states(states_sign):
    full_states = [1, 4 / 3, 1.5, 2, 2.5, 3, 4, 5, 6]
    full_labels = ['1', '4/3', '3/2', '2', '5/2', '3', '4', '5', '6']
    full_colors = (
        '#56B4E9',
        '#0072B2',
        '#009E73',
        '#E69F00',
        '#F0E442',
        '#D55E00',
        '#999999',
        '#505050',
        '#CC79A7')
    if states_sign == 'int_5':
        states = [0, 3, 5, 6, 7]
    elif states_sign == 'int_6':
        states = [0, 3, 5, 6, 7, 8]
    elif states_sign == 'full_5':
        states = [0, 2, 3, 5, 6, 7]
    elif states_sign == 'full_5_and_6':
        states = [0, 2, 3, 5, 6, 7, 8]
    elif states_sign == 'full_6_but_1.33':
        states = [0, 2, 3, 4, 5, 6, 7, 8]
    elif states_sign == 'full_6_but_2.5':
        states = [0, 1, 2, 3, 5, 6, 7, 8]
    elif states_sign == 'full_6':
        states = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    else:
        states = []
    return ([full_states[x] for x in states],
            [full_labels[x] for x in states],
            [full_colors[x] for x in states])

sns.set(font_scale=1.4, style="ticks", font="lato",
        # palette=('#56B4E9', '#009E73', '#F0E442'))
        # palette=('#7570b3', '#d95f02', '#1b9e77'))
        palette=('#56B4E9', '#E69F00', '#009E73', '#D55E00', '#CC79A7', '#F0E442'))
# palette=('#1f77b4', '#2ca02c', '#ff7f0e'))
sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})
plt.rcParams['font.weight'] = "medium"
plt.rcParams['axes.labelweight'] = 'medium'
plt.rcParams['figure.titleweight'] = 'medium'
plt.rcParams['axes.titleweight'] = 'medium'
plt.rcParams['figure.figsize'] = 8, 10*2/3
plt.rcParams["legend.framealpha"] = 0.5
plt.rcParams['axes.xmargin'] = 0
plt.rcParams['axes.ymargin'] = 0
plt.rcParams["legend.framealpha"] = 1

cells = ('ALL', 'K562', 'MCF7', 'A549', 'HCT116', '22RV1', 'Other')
markers = ('s', '*', 'v', '1', 'o')
state_ss = ['int_6', 'full_5_and_6', 'full_6_but_1.33', 'full_6_but_2.5', 'full_6']
short_s = dict(zip(state_ss, ['int_6', 'full_5\nand_6', 'full_6\nbut_1.33', 'full_6\nbut_2.5', 'full_6']))
state_ss = ['full_5_and_6']
short_s['COSMIC'] = 'COSMIC'
# colors = ('#E69F00', '#009E73', '#D55E00', '#CC79A7', '#F0E442', '#56B4E9')
all_states, all_labels, all_colors = get_states('full_6')
for cell_sign in ['K562']:#cells:
    print(cell_sign)
    fig, (*axs,) = plt.subplots(2, 1)
    # axs = [ax for tup in axs for ax in tup]
    for state_s, ax in zip(['COSMIC'] + state_ss, axs):
        # model = 'CAIC@{}@{}'.format(state_s if state_s != 'COSMIC' else 'full_6', 4)
        model='CAIC'
        print(model)
        states, labels, colors = get_states(state_s if state_s != 'COSMIC' else 'full_6')
        # sns.set_palette(colors)
        t = {}
        min_tr = {}
        max_tr = {}
        for BAD in states:
            if cell_sign == 'ALL':
                t[BAD] = pd.read_table(os.path.expanduser(
                    'D:\Sashok\Desktop/counts/counts_deltaqm_{}_{:.2f}.tsv'.format(model, BAD)))
            else:
                t[BAD] = pd.read_table(os.path.expanduser(
                    'D:\Sashok\Desktop/counts/counts_deltaqm_{}_{}_{:.2f}.tsv'.format(cell_sign, model, BAD)))
            t[BAD].replace(1.3333333333333337, 4 / 3, inplace=True)
            min_tr[BAD] = t[BAD]['threshold'].min()
            max_tr[BAD] = t[BAD]['threshold'].max()
        cosmic_states = list(t[1]['COSMIC'].unique())
        SNPS_DIST = {BAD: 0 for BAD in list(set(all_states) | set(cosmic_states))}
        if state_s == 'COSMIC':
            print(cosmic_states)
            for BAD in cosmic_states:
                SNPS_DIST[BAD] = t[1][(t[1]['COSMIC'] == BAD) & (t[1]['threshold'] == min_tr[1])]['counts'].sum()
            ALL = t[1][t[1]['threshold'] == min_tr[1]]['counts'].sum()

            #!!!
            print(SNPS_DIST[6])
            SNPS_DIST[6] = t[1][(t[1]['COSMIC'] >= 6) & (t[1]['threshold'] == min_tr[1])]['counts'].sum()
            print(SNPS_DIST[6])
            for x in SNPS_DIST.keys():
                if x > 6:
                    SNPS_DIST[x] = 0

            states_to_draw = sorted(list(SNPS_DIST.keys()))

            sum_vals = sum(SNPS_DIST.values())
            if ALL != sum_vals:
                print('asd', ALL, sum_vals)
            snps_list = [SNPS_DIST[x] / ALL for x in states_to_draw]

            states_to_draw, snps_list = \
                zip(*[(state, y) for state, y in zip(states_to_draw, snps_list) if y >= 0.005 or state in all_states])
        else:
            actual_seg_tr = dict(
                zip(states, [min((x for x in list(t[BAD]['threshold'].unique()) if x >= 0), default=0) for BAD in states]))
            for BAD in states:
                SNPS_DIST[BAD] = t[BAD][t[BAD]['threshold'] == actual_seg_tr[BAD]]['counts'].sum()

            ALL = sum(SNPS_DIST[BAD] for BAD in states)

            sum_vals = sum(SNPS_DIST.values())
            if ALL != sum_vals:
                print('asd', ALL, sum_vals)
            snps_list = [SNPS_DIST[x] / ALL for x in states_to_draw]

        print(dict(zip(states_to_draw, snps_list)), sum(snps_list))

        all_colors_dict = {state: color for state, color in zip(all_states, all_colors)}
        all_labels_dict = {state: label for state, label in zip(all_states, all_labels)}

        colors_to_draw = [all_colors_dict[s] if s in all_colors_dict else '#FFFFFF' for s in states_to_draw]
        labels_to_draw = [all_labels_dict[s] if s in all_labels_dict else get_state_label('{:.2f}'.format(s)) for s in states_to_draw]

        ax.bar(x=list(range(len(states_to_draw))), height=snps_list,
               color=colors_to_draw, edgecolor='black', tick_label=labels_to_draw)
        ax.set_ylim(0, 1)
        ax.set_ylabel(short_s[state_s])
        ax.grid()
        # if state_s == 'COSMIC':
        #     ax.set_title('BAD distribution on SNPs\n{}'.format(cell_sign))
    plt.suptitle('BAD distribution on SNPs\n{}'.format(cell_sign))
    print('D:\Sashok\Desktop/PR_SUM/Dist_{}.png'.format(cell_sign))
    plt.savefig(os.path.expanduser('D:\Sashok\Desktop/PR_SUM/Dist_{}.png'.format(cell_sign)), dpi=300)
