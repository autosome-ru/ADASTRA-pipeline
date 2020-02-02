import sys
import os
from matplotlib import pyplot as plt

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.helpers import ChromPos

fig = plt.figure(figsize=(24, 8))


class ChromPos:
    chrs = {'chrX', 'chrY'}
    for i in range(1, 22):
        chrs.add('chr' + str(i))


def plot_ploidy(file, C, i, cnv=False, name=None):
    print(file)
    if cnv:
        for line in file:
            if line[0] == '#':
                continue
            line = line.strip().split(',')
            # if int(line[4]) in {4,6,8} or line[3] == '0': continue
            if line[0] != name:
                continue
            if 'chr' + line[4] not in ChromPos.chrs:
                continue
            if 'chr' + line[4] != chr:
                continue
            if int(line[10]) == 0:
                y = 0
            else:
                y = int(line[11]) / int(line[10]) - 1

            try:
                start = int(line[5])
                end = int(line[6])
            except ValueError:
                continue
            if y >= 1:
                plt.hlines(y=y, xmin=start, xmax=end, colors=C, linewidth=7)
    else:
        for line in file:
            if line[0] == '#':
                continue
            line = line.split()
            if line[0] != chr:
                continue
            start = int(line[1])
            end = int(line[2])
            y = float(line[3])
            if y != 0:
                plt.hlines(y=y, xmin=start, xmax=end, colors=C, linewidth=6 - 2 * i)
    # plt.vlines(x=end, ymin=0, ymax=10, linestyle='--', linewidth=1, color=C)


# name = 'HCT-116_colon_carcinoma!_labs_richard-myers___biosamples_ENCBS389ENC_'
name = 'HCT-116_colon_carcinoma!_labs_michael-snyder___biosamples_ENCBS626JHZ_'
snps_name = os.path.expanduser('~/Documents/ASB/simulation/' + name + '.tsv')
ploidy_name = os.path.expanduser('~/Documents/ASB/simulation/' + name + '_ploidy.tsv')
cnv_line = 'HCT-116'


for chr in ChromPos.chrs:
    print(chr)
    fig, ax = plt.subplots()
    with open('/home/sashok/Documents/ASB/Cell_lines/cell_lines_copy_number.csv', 'r') as file:
        plot_ploidy(file, 'C' + str((0) % 10), 0, cnv=True, name=cnv_line)

    with open(ploidy_name, 'r') as file:
        plot_ploidy(file, 'C' + str((1) % 10), 1)

    x = []
    y = []
    colors = []
    with open(snps_name, 'r') as snps:
        for line in snps:
            if line[0] == '#':
                continue
            line = line.split()
            if line[0] != chr:
                continue
            if line[2] == '.':
                continue
            pos = int(line[1])
            ref = int(line[5])
            alt = int(line[6])
            cov = ref + alt
            cm = 100
            x.append(pos)
            # colors.append([0.2+0.8*(1-cov/cm), 0.2, 0.2+0.8*cov/cm])
            # colors.append(cov)
            if cov < 8:
                colors.append('Grey')
            elif cov < 30:
                colors.append('Grey')
            elif cov < 50:
                colors.append('Green')
            else:
                colors.append('Purple')
            if min(ref, alt) != 0:
                y.append(max(ref, alt) / min(ref, alt))
            else:
                y.append(0)
    if not x:
        continue
    print(len(x), len(y), min(y), max(y))

    plt.scatter(x, y, c=colors, s=1)
    plt.plot([1, 1], [1, 1], color='C0', label='CNV')
    plt.plot([1, 1], [1, 1], color='C1', label='estimations')
    plt.legend()
    plt.grid(True)
    plt.xlabel('position, bp')
    plt.ylabel('1/MAF-1')
    plt.title(chr + '   ' + 'Corrected-1.5-CAIC-NEW_CG')
    plt.savefig(os.path.expanduser('~/CPPloidy/' + chr + '_' + name + '.png'))
    plt.close(fig)
