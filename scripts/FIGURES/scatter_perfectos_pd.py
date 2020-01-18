import sys
import numpy as np
import os
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


def get_color(row):
    if abs(row['log_fc']) < np.log10(fc_tr) or \
            (abs(row['fdrp_by_ref']) > fdr_tr and abs(row['fdrp_by_alt']) > fdr_tr)\
            or (row['motif_log_pref'] < -np.log10(perf_tr)) or (row['motif_log_palt'] < -np.log10(perf_tr)):
        return 'grey'
    if row['log_fc'] * row['log_pv'] > 0:
        return 'blue'
    else:
        return 'red'


if __name__ == '__main__':
    filepath = os.path.expanduser("~/CTCF_HUMAN_fc.tsv")
    pt = pd.read_table(filepath)

    perf_tr = 0.00005
    fc_tr = 2
    fdr_tr = 0.05
    fix = ''

    pt = pt[~(pt['fdrp_by_alt'].isnull() | pt['fdrp_by_ref'].isnull())]
    # pt = pt[(pt['fdrp_by_alt'] <= fdr_tr) | (pt['fdrp_by_ref'] <= fdr_tr)]
    print(pt)
    pt['log_pv'] = (np.log10(
        pt['fdrp_by_ref']) - np.log10(pt['fdrp_by_alt']))
    pt['log_fc'] = pt['fold_change']
    pt['col'] = pt.apply(lambda x: get_color(x), axis=1)

    print(pt.info())

    blue = len(pt[pt.col == 'blue'].index)
    red = len(pt[pt.col == 'red'].index)
    grey = len(pt[pt.col == 'grey'].index)

    print(pt['log_pv'])

    fig, ax = plt.subplots(figsize=(10, 8))
    plt.scatter(x=list(pt['log_pv']), y=list(pt['log_fc']), c=pt['col'], s=5)
    plt.grid(True)
    plt.title(
        '\n'.join(['perfectos_pv <= {}', 'fold_change >= {}', 'fdr <= {}', fix]).format(perf_tr, fc_tr, fdr_tr))
    plt.xlabel('signed -log10 fdr_p')
    plt.ylabel('log10 prefectos foldchange')
    label = 'blue/red: {0}/{1}({2}%),\ngrey/all={3}%'.format(blue, red,
                                                             round(blue / (blue + red) * 100, 1),
                                                             round(grey / (grey + blue + red) * 100, 1))
    plt.text(x=max(pt.log_pv) / 5, y=min(pt.log_fc) / 2, s=label)
    plt.text(x=min(pt.log_pv) * 4 / 5, y=min(pt.log_fc) / 2,
             s=' log P-value treshold: {},\nlog FC treshold: {}'.format(round(fdr_tr, 1), round(fc_tr, 1)))

    plt.savefig(os.path.expanduser("~/{}_p_tr={:.2f}_fc_tr={:.2f}_fdr_tr={:.2f}.png".format(filepath.split('/')[-1],
                                                                                perf_tr, fc_tr, fdr_tr)))
