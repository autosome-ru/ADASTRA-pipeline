import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, ticker
import seaborn as sns

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.helpers import states

sns.set(font_scale=1.4, style="ticks", font="lato", palette=(
"#f15854", "#faa43a", "#e5d00d", "#60bd68", "#5da5da", "#f17cb0", "#975597", "#b2912f", "#aaaaaa", "#4d4d4d"))
# sns.set(palette=('#fcfbfd', '#efedf5', '#dadaeb', '#bcbddc', '#9e9ac8', '#807dba', '#6a51a3', '#54278f', '#3f007d'))
# sns.set(palette=(
# "#f15854", "#faa43a", "#e5d00d", "#60bd68", "#5da5da", "#f17cb0", "#975597", "#b2912f", "#aaaaaa", "#4d4d4d"))
sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})
plt.rcParams['font.weight'] = "medium"
plt.rcParams['axes.labelweight'] = 'medium'
plt.rcParams['figure.titleweight'] = 'medium'
plt.rcParams['axes.titleweight'] = 'medium'
plt.rcParams['figure.figsize'] = 6, 5
plt.rcParams["legend.framealpha"] = 1
plt.rcParams['axes.xmargin'] = 0
plt.rcParams['axes.ymargin'] = 0
plt.rcParams["legend.framealpha"] = 0.5

t = pd.read_table(os.path.expanduser('~/counts_qm.tsv'))

t.replace(1.3333333333333337, 4/3, inplace=True)
# for x in t['COSMIC'].unique():
#     print(x, t[t['COSMIC'] == x]['counts'].sum())

x = [x for x in range(101)] + [x for x in range(110, 201, 10)]

TP = {}
FP = {}
Precision = {}  # TP / (TP + FP)
Recall = {}  # TPR = TP / P
FPR = {}  # FP / N
for BAD in states:
    TP[BAD] = {}
    FP[BAD] = {}
    Precision[BAD] = {}
    Recall[BAD] = {}
    FPR[BAD] = {}

All = {}
Global_precision = {}

t = t[t['BAD'].isin(states) & t['COSMIC'].isin(states)]
l = len(states)

P = {}
N = {}
for BAD in states:
    P[BAD] = t[t['COSMIC'] == BAD]['counts'].sum()
    N[BAD] = t[t['COSMIC'] != BAD]['counts'].sum()

print('starto')

for tr in x:
    s = t[t['threshold'] == tr].copy()
    All[tr] = s['counts'].sum()
    Global_precision[tr] = s[s['BAD'] == s['COSMIC']]['counts'].sum() / All[tr]
    for BAD in states:
        # print(tr, BAD)
        TP[BAD][tr] = s[(s['COSMIC'] == BAD) & (s['BAD'] == BAD)]['counts'].sum()
        FP[BAD][tr] = s[(s['COSMIC'] != BAD) & (s['BAD'] == BAD)]['counts'].sum()
        Precision[BAD][tr] = TP[BAD][tr] / (TP[BAD][tr] + FP[BAD][tr])
        Recall[BAD][tr] = TP[BAD][tr] / P[BAD]
        FPR[BAD][tr] = FP[BAD][tr] / N[BAD]

# s[((s['COSMIC'] == BAD) & (s['BAD'] == BAD)) | ((s['COSMIC'] != BAD) & (s['BAD'] != BAD))]['counts'].sum() / All)

print('starto2')

# P[tr] R[tr]
for metric, label in (Precision, 'Precision'), (Recall, 'Recall'):
    fig, ax = plt.subplots()
    for BAD in states:
        plt.plot(x, [metric[BAD][tr] for tr in x], label='{:.2f}'.format(BAD))
    plt.plot(x, [Global_precision[tr] for tr in x], label='All', color='black')
    # ax.axhline(y=((l - 1) ** 2 + 1) / ((l - 1) ** 2 + 1 + 2 * (l - 1)), color='black', linestyle='--')
    ax.axhline(y=1 / l, color='black', linestyle='--')

    plt.legend(loc='lower right')
    plt.grid(True)
    plt.xlabel('SNP qual >= x')
    plt.ylabel(label)

    plt.savefig(os.path.expanduser('~/AC_5/figa_{}.png'.format(label)))
    plt.close(fig)

print('starto3')

# PR-curve
fig, ax = plt.subplots()
for BAD in states:
    ax.plot([Recall[BAD][tr] for tr in x], [Precision[BAD][tr] for tr in x], label='{:.2f}'.format(BAD))
ax.axhline(y=1 / l, color='black', linestyle='--')
# ax.axvline(x=1 / l, color='black', linestyle='--')
plt.legend(loc='lower center')
plt.grid(True)
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.savefig(os.path.expanduser('~/AC_5/figa_PR.png'))
plt.close(fig)

print('starto4')

# ROC
fig, ax = plt.subplots()
for BAD in states:
    ax.plot([Recall[BAD][tr] for tr in x], [FPR[BAD][tr] for tr in x], label='{:.2f}'.format(BAD))
ax.axhline(y=1 / l, color='black', linestyle='--')
ax.axvline(x=1 / l, color='black', linestyle='--')
plt.legend(loc='lower right')
plt.grid(True)
plt.xlabel('TPR')
plt.ylabel('FPR')
plt.savefig(os.path.expanduser('~/AC_5/figa_ROC.png'))

plt.show()
plt.close(fig)
