import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, ticker
import seaborn as sns

sys.path.insert(1, "/home/abramov/ASB-Project")


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


sns.set(font_scale=1.4, style="ticks", font="lato")

markers = ('H', '*', 'v', '1', 'o')
state_ss = ('int_6', 'full_5_and_6', 'full_6_but_1.33', 'full_6_but_2.5', 'full_6')

df_cell = pd.read_table(os.path.expanduser('D:\Sashok/Desktop/AC_5/cell_wise_stats.tsv'))
df_all = pd.read_table(os.path.expanduser('D:\Sashok/Desktop/AC_5/all_stats.tsv'))

all_states, all_labels, _ = get_states('full_6')

for cell_sign in ('ALL', 'K562', 'MCF7', 'A549', 'HCT116', '22RV1', 'Other'):
    if cell_sign == 'ALL':
        df = df_all
    else:
        df = df_cell
    fig, ax = plt.subplots()
    ax.grid(True)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_title('PR summary\n{}'.format(cell_sign))
    for state_sign, marker in zip(state_ss, markers):
        states, labels, colors = get_states(state_sign)
        for mult in range(4, 5):
            for BAD, label, color in zip(states, labels, colors):
                precision = df[(df['states'] == state_sign) & (df['multiplier'] == mult)
                    ]['Precision@{}@BAD={}'.format(cell_sign, label)]
                recall = df[(df['states'] == state_sign) & (df['multiplier'] == mult)
                    ]['Recall@{}@BAD={}'.format(cell_sign, label)]
                ax.scatter(x=[recall], y=[precision], marker=marker, color=color, s=10*mult, alpha=0.5)
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    legend_elements = [plt.scatter()]
    legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='BAD', frameon=True, fancybox=False,
                       labels=all_labels)
    legend.get_frame().set_edgecolor('black')
    plt.savefig(os.path.expanduser('D:\Sashok\Desktop/PR_SUM/PR_SUM_{}.png'.format(cell_sign)),
                dpi=300)
    plt.close(fig)
