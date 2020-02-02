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

t = pd.read_table(os.path.expanduser('~/counts_m.tsv'))

const = 'COSMIC'
t.replace(1.3333333333333337, 4/3, inplace=True)
# for x in t['COSMIC'].unique():
#     print(x, t[t['COSMIC'] == x]['counts'].sum())

x = []
y = {}
for const in 'BAD', 'COSMIC':
    y[const] = {}
    for BAD in states:
        y[const][BAD] = []
y['overall'] = []

t = t[t['BAD'].isin(states)]
l = len(states)

for tr in range(10, 101):
    x.append(tr)
    s = t[t['threshold'] == tr].copy()
    All = s['counts'].sum()
    y['overall'].append(s[s['BAD'] == s['COSMIC']]['counts'].sum() / All)
    for const in 'BAD', 'COSMIC':
        for BAD in states:
            All_cur = s[(s[const] == BAD)]['counts'].sum()
            if All_cur == 0:
                y[const][BAD].append(0)
            else:
                y[const][BAD].append(s[(s['COSMIC'] == BAD) & (s['BAD'] == BAD)]['counts'].sum() / All_cur)

# s[((s['COSMIC'] == BAD) & (s['BAD'] == BAD)) | ((s['COSMIC'] != BAD) & (s['BAD'] != BAD))]['counts'].sum() / All)

for const in 'BAD', 'COSMIC':
    fig, ax = plt.subplots()
    for BAD in states:
        plt.plot(x, y[const][BAD], label='{:.2f}'.format(BAD))
    plt.plot(x, y['overall'], label='All', color='black')
    # ax.axhline(y=((l - 1) ** 2 + 1) / ((l - 1) ** 2 + 1 + 2 * (l - 1)), color='black', linestyle='--')
    ax.axhline(y=1 / l, color='black', linestyle='--')

    # ax.set_facecolor('white')
    plt.legend(loc='lower right')
    plt.grid(True)
    plt.xlabel('SNP coverage >= x')
    plt.ylabel('{} (const {})'.format({'COSMIC': 'Recall', 'BAD': 'Precision'}[const], const))

    plt.savefig(os.path.expanduser('~/AC_5/figa_{}.png'.format({'COSMIC': 'Recall', 'BAD': 'Precision'}[const])))
    plt.close(fig)

fig, ax = plt.subplots()
for BAD in states:
    ax.plot(y['COSMIC'][BAD], y['BAD'][BAD], label='{:.2f}'.format(BAD))
ax.axhline(y=1 / l, color='black', linestyle='--')
ax.axvline(x=1 / l, color='black', linestyle='--')
plt.legend(loc='lower right')
plt.grid(True)
plt.xlabel('Recall (const COSMIC)')
plt.ylabel('Precision (const BAD)')
plt.savefig(os.path.expanduser('~/AC_5/figa_PR.png'))

plt.show()
plt.close(fig)