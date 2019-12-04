import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


def get_color(pv, fc):
    if abs(fc) <= max(fc_tr, 0) or pv <= max(fdr_tr, 0):
        return 'grey'
    if fc * pv > 0:
        return 'blue'
    else:
        return 'green'


pt = pd.read_table(sys.argv[1])

perf_tr = 0.05
fc_tr = np.log10(2)
fdr_tr = -np.log10(0.05)
fix = ''

pt = pt[(pt.prefectos_p1 <= perf_tr) & (pt.perfectos_p2 <= perf_tr)]
pt.log_fc = np.log10(pt.perfectos_fc)
pt.log_pv = np.sign(pt['fdr_alt' + fix] - pt['fdr_ref' + fix]) * np.log10(min(pt['fdr_alt' + fix], pt['fdr_ref' + fix]))
pt['col'] = pt.apply(lambda x: get_color(x.log_pv, x.log_fc), axis=1)

blue = len(pt[pt.col == get_color(1000, 1000)].index)
red = len(pt[pt.col == get_color(1000, -1000)].index)
grey = len(pt[pt.col == get_color(0, 0)].index)

fig, ax = plt.subplots(figsize=(10, 8))
sns.scatterplot(data=pt, x='log_pv', y='log_fc', color='col')
plt.grid(True)
plt.title(
    '\n'.join(['perfectos_pv <= {}', 'fold_change >= {}', 'fdr <= {}', fix]).format(perf_tr, 10 ** fc_tr, 10 ** fdr_tr))

label = 'blue/red: {0}/{1}({2}%),\ngrey/all={3}%'.format(blue, red,
                                                         round(blue / (blue + red) * 100, 1),
                                                         round(grey / (grey + blue + red) * 100, 1))
plt.text(x=max(pt.log_pv) / 5, y=min(pt.log_fc) / 2, s=label)
plt.text(x=min(pt.log_pv) * 4 / 5, y=min(pt.log_fc) / 2,
         s='P-value treshold: {},\nFC treshold: {}'.format(round(fdr_tr, 1), round(fc_tr, 1)))

plt.savefig('_'.join([sys.argv[1].split('/')[-1], perf_tr, fc_tr, fdr_tr, fix]) + '.png')
