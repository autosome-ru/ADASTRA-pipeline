import sys
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats


def unpack(line, header):
    line = line.strip().split('\t')
    assert len(line) == len(header)
    data = dict(zip(header, line))

    for x in ['sBAD',  'm_fpref', 'm_fpalt',
              'm_stpref', 'm_stpalt', 'perfectos_p1', 'perfectos_p2', 'perfectos_fc']:
        if x in data:
            data[x] = float(data[x])

    for x in ['m_datasets', 'motif_pos']:
        if x in data:
            if data[x] != '.':
                data[x] = int(data[x])

    return data


path = sys.argv[1]
with open(path, 'r') as file:
    pv = []
    fc = []
    colors = []
    Ivans_nightmare = 0
    grey = 0
    blue = 0
    red = 0
    perf_tr = 0.05
    # header по умолчанию
    header = 'chr    pos     ID      ref     alt     m_callers       m_ploidy        m_q     m_dipq  m_segc  m_datasets      m_hpref m_hpalt m_fpref m_fpalt m_stpref        m_stpalt        perfectos_p1    perfectos_p2    perfectos_fc    motif_pos       orientation'.split(
        '\t')
    for line in file:
        if line == "":
            continue
        if line[0] == '#':
            header = line[1:].strip().split('\t')
            continue
        args = unpack(line, header)
        if args['perfectos_p1'] >= perf_tr and args['perfectos_p2'] >= perf_tr:
            continue

        if args['fdrp_ref'] < 1:
            if args['fdrp_ref'] == 0:
                continue
            if not args['fdrp_alt'] == 1:
                Ivans_nightmare += 1
                continue
            pv.append(np.log10(args['fdrp_ref']))
            fc.append(np.log10(args['perfectos_fc']))
            continue
        if args['fdrp_alt'] < 1:
            if args['fdrp_alt'] == 0:
                continue
            if not args['fdrp_alt'] == 1:
                Ivans_nightmare += 1
                continue
            pv.append(-np.log10(args['fdrp_alt']))
            fc.append(np.log10(args['perfectos_fc']))
            continue

    mp = max([abs(x) for x in pv])
    fc_tr = np.log10(2)  # perfectos treshold
    pv_tr = 1

    for i, p in enumerate(pv):
        f = fc[i]
        if abs(f) < fc_tr or abs(p) < pv_tr:
            colors.append('grey')
            grey += 1
        elif f * p > 0:
            colors.append('blue')
            blue += 1
        else:
            colors.append('red')
            red += 1
    filt = "path"
    label = 'blue/red: {0}/{1}({2}%),\ngrey/all={3}%'.format(blue, red,
                                                             round(blue / (blue + red) * 100, 1),
                                                             round(grey / (grey + blue + red) * 100, 1))
    print(stats.spearmanr(pv, fc))
    print(stats.pearsonr(pv, fc))
    print(label)
    title = path.split("/")[-1].split(".")[0].split("_")[0] + " " + path.split("/")[-1].split(".")[0].split("_")[3]
    print(title)
    plt.scatter(pv, fc, s=5, c=colors)
    plt.title(title)
    plt.text(x=max(pv) / 5, y=min(fc) / 2, s=label)
    plt.text(x=min(pv) * 4 / 5, y=min(fc) / 2,
             s='P-value treshold: {},\nFC treshold: {}'.format(round(pv_tr, 1), round(fc_tr, 1)))
    plt.xlabel('signed -1*log10(fdr))')
    plt.ylabel('perfectos fc')
    plt.grid(True)
    plt.show()
