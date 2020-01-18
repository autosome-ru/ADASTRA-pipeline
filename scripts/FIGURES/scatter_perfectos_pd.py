import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


def get_color(pv, fc):
    if abs(fc) < fc_tr or abs(pv) < fdr_tr:
        return 'grey'
    if fc * pv > 0:
        return 'blue'
    else:
        return 'red'


if __name__ == '__main__':
    pt = pd.read_table(sys.argv[1])

    perf_tr = 1.3
    fc_tr = np.log10(2)
    fdr_tr = -np.log10(0.05)
    fix = ''

    pt = pt[(pt['motif_log_pref'] >= perf_tr) & (pt['motif_log_palt'] >= perf_tr)]
    pt = pt[(0 < pt['fdrp_by_ref' + fix]) & (0 < pt['fdrp_by_alt' + fix])]

    pt['log_fc'] = pt['fold_change']
    pt['log_pv'] = np.sign(pt['fdrp_by_alt' + fix] - pt['fdrp_by_ref' + fix]) * (np.log10(
        pt['fdrp_by_alt']) - np.log10(pt['fdrp_by_ref']))
    pt['col'] = pt.apply(lambda x: get_color(x['log_pv'],
                                             x['log_fc']), axis=1)

    print(pt.info())

    blue = len(pt[pt.col == get_color(1000, 1000)].index)
    red = len(pt[pt.col == get_color(1000, -1000)].index)
    grey = len(pt[pt.col == get_color(0, 0)].index)

    print(pt['log_pv'])

    fig, ax = plt.subplots(figsize=(10, 8))
    plt.scatter(x=list(pt['log_pv']), y=list(pt['log_fc']), c=pt['col'], s=5)
    plt.grid(True)
    plt.title(
        '\n'.join(['perfectos_pv <= {}', 'fold_change >= {}', 'fdr <= {}', fix]).format(perf_tr, 10 ** fc_tr, round(10 ** (-fdr_tr), 2)))
    plt.xlabel('signed -log10 fdr_p')
    plt.ylabel('log10 prefectos foldchange')
    label = 'blue/red: {0}/{1}({2}%),\ngrey/all={3}%'.format(blue, red,
                                                             round(blue / (blue + red) * 100, 1),
                                                             round(grey / (grey + blue + red) * 100, 1))
    plt.text(x=max(pt.log_pv) / 5, y=min(pt.log_fc) / 2, s=label)
    plt.text(x=min(pt.log_pv) * 4 / 5, y=min(pt.log_fc) / 2,
             s=' log P-value treshold: {},\nlog FC treshold: {}'.format(round(fdr_tr, 1), round(fc_tr, 1)))

    plt.savefig('_'.join(map(str, [sys.argv[1].split('/')[-1], perf_tr, fc_tr, fdr_tr, fix])) + '.png')
