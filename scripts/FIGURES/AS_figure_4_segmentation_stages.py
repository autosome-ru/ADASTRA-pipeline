import sys
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt, ticker
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.helpers import ChromPos

sns.set(font_scale=1.2, style="ticks", font="lato", palette=('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2',
                                                             '#D55E00', '#CC79A7'))
# sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})
plt.rcParams['font.weight'] = "medium"
plt.rcParams['axes.labelweight'] = 'medium'
plt.rcParams['figure.titleweight'] = 'medium'
plt.rcParams['axes.titleweight'] = 'medium'
plt.rcParams['figure.figsize'] = 14, 6
plt.rcParams["legend.framealpha"] = 1
plt.rcParams['axes.xmargin'] = 0
plt.rcParams['axes.ymargin'] = 0
plt.rcParams["legend.framealpha"] = 1

# name = 'HCT-116_colon_carcinoma!_labs_richard-myers___biosamples_ENCBS389ENC_'
# name = 'HCT-116_colon_carcinoma!_labs_michael-snyder___biosamples_ENCBS626JHZ_'
# snps_name = os.path.expanduser('~/Documents/ASB/simulation/' + name + '.tsv')
# ploidy_name = os.path.expanduser('~/Documents/ASB/simulation/' + name + '_ploidy.tsv')

# snps_name = os.path.expanduser('~/cherry_BAD/K562__myelogenous_leukemia_!_labs_xiang-dong-fu___biosamples_ENCBS074NGX_.tsv')
# ploidy_name = os.path.expanduser('~/K562_BAD_Segments/K562_myelogenous_leukemia!_labs_xiang-dong-fu___biosamples_ENCBS074NGX__ploidy.tsv')

snps_name = os.path.expanduser('~/DataForFigures/K562__myelogenous_leukemia_!_labs_michael-snyder___biosamples_ENCBS725WFV_.tsv')
ploidy_name = os.path.expanduser('~/DataForFigures/K562__myelogenous_leukemia_!_labs_michael-snyder___biosamples_ENCBS725WFV__ploidy.tsv')


cosmic_name = os.path.expanduser('~/DataForFigures/cell_lines_copy_number.csv')
cnv_line = 'K-562'

snps = pd.read_table(snps_name, header=None)
snps.columns = ['chr', 'pos', 'ID', 'ref', 'alt', 'ref_c', 'alt_c']
snps = snps[snps['ID'] != '.']
snps['AD'] = snps[['ref_c', 'alt_c']].max(axis=1) / snps[['ref_c', 'alt_c']].min(axis=1)
snps['cov'] = snps['ref_c'] + snps['alt_c']
snps['log_cov'] = np.log10(snps['cov'])

ploidy = pd.read_table(ploidy_name)

# chrs = ('chr10', 'chr17')
chrs = ('chr12', 'chr6')
BAD_color = '#0072B266'
BAD_color_1 = '#0072B2CC'
COSMIC_color = '#D55E00'
BAD_lw = 10
COSMIC_lw = 4
y_min = 0.8
y_max = 6
delta_y = 0.05

# First step
fig, (*axs,) = plt.subplots(len(chrs), 1)
fig.tight_layout(pad=1.5)
plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))

for chr, ax in zip(chrs, axs):
    chr_ploidy = ploidy[ploidy['#chr'] == chr]
    chr_snps = snps[snps['chr'] == chr].copy()
    chr_snps['AD'] = chr_snps['AD'].apply(lambda y: y_max - delta_y if y > y_max else y)

    bar_positions = []
    bar_widths = []
    bar_colors = []
    vd = 1 / 500

    borders_to_draw = []
    last_end = 1
    for index, (pl_chr, start, end, BAD, *values) in chr_ploidy.iterrows():
        if start != last_end:
            if last_end == 1:
                borders_to_draw += [start - ChromPos.chrs[chr] * vd]
            else:
                borders_to_draw += [last_end + ChromPos.chrs[chr] * vd, start - ChromPos.chrs[chr] * vd]
            bar_colors.append('#AAAAAA')
        last_end = end
        bar_colors.append('C2')
    if last_end != ChromPos.chrs[chr] + 1:
        borders_to_draw += [last_end + ChromPos.chrs[chr] * vd]
        bar_colors.append('#AAAAAA')

    reduced_bar_colors = []
    for i, color in enumerate(bar_colors):
        if i == 0 or bar_colors[i - 1] != color:
            reduced_bar_colors.append(color)

    borders_for_bars = [1] + borders_to_draw + [ChromPos.chrs[chr] + 1]
    for i in range(len(borders_for_bars) - 1):
        bar_positions.append((borders_for_bars[i] + borders_for_bars[i + 1]) / 2)
        bar_widths.append(borders_for_bars[i + 1] - borders_for_bars[i])

    for border in borders_to_draw:
        ax.axvline(x=border, ymin=0, ymax=0.5, linestyle='--', color='C2')

    ax.scatter(x='pos', y='AD', color='#555555', data=chr_snps, alpha=0.5, s=2, vmin=10, vmax=30)
    ax.set_xlim(0, ChromPos.chrs[chr])
    ax.set_ylim(y_min, y_max)
    ax.grid(which='major', axis='both')
    ax.set_xticklabels([])
    ax.set_yticks(list(range(1, int(y_max) + 1)))
    ax.text(0.99, 0.95, 'K562 {}'.format(chr),
            horizontalalignment='right',
            verticalalignment='top',
            transform=ax.transAxes)
    ax.set_ylabel('AD')
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="10%", pad=0.05)
    cax.get_yaxis().set_ticks([])
    cax.set_xlim(1, ChromPos.chrs[chr])
    cax.set_ylim(0, 1)
    cax.bar(bar_positions, [1] * len(bar_positions), bar_widths, color=reduced_bar_colors, linewidth=0)

cax.set_xlabel('Chromosome position, bp')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)

plt.savefig(os.path.expanduser('~/AC_4/Figure_AS_4_stage1.svg'), dpi=300)
plt.savefig(os.path.expanduser('~/AC_4/Figure_AS_4_stage1.png'), dpi=300)
# plt.show()
plt.close(fig)

# BAD step
fig, (*axs,) = plt.subplots(len(chrs), 1)
fig.tight_layout(pad=1.5)
plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))

for chr, ax in zip(chrs, axs):
    chr_ploidy = ploidy[ploidy['#chr'] == chr]
    chr_snps = snps[snps['chr'] == chr].copy()
    chr_snps['AD'] = chr_snps['AD'].apply(lambda y: y_max - delta_y if y > y_max else y)

    bar_positions = []
    bar_widths = []
    bar_colors = []
    vd = 1 / 500

    BADs = []

    borders_to_draw = []
    segmentation_borders = []
    last_end = 1
    for index, (pl_chr, start, end, BAD, *values) in chr_ploidy.iterrows():
        if start != last_end:
            if last_end == 1:
                borders_to_draw += [start - ChromPos.chrs[chr] * vd]
                segmentation_borders += [start - ChromPos.chrs[chr] * vd]
            else:
                borders_to_draw += [last_end + ChromPos.chrs[chr] * vd, start - ChromPos.chrs[chr] * vd]
                segmentation_borders += [last_end + ChromPos.chrs[chr] * vd, start - ChromPos.chrs[chr] * vd]
            bar_colors.append('#AAAAAA')
            BADs.append(None)
        else:
            if last_end != 1:
                segmentation_borders += [last_end]
        last_end = end
        bar_colors.append('C2')
        BADs.append(BAD)
    if last_end != ChromPos.chrs[chr] + 1:
        borders_to_draw += [last_end + ChromPos.chrs[chr] * vd]
        segmentation_borders += [last_end + ChromPos.chrs[chr] * vd]
        bar_colors.append('#AAAAAA')
        BADs.append(None)

    reduced_bar_colors = []
    for i, color in enumerate(bar_colors):
        if i == 0 or bar_colors[i - 1] != color:
            reduced_bar_colors.append(color)

    borders_for_bars = [1] + borders_to_draw + [ChromPos.chrs[chr] + 1]
    for i in range(len(borders_for_bars) - 1):
        bar_positions.append((borders_for_bars[i] + borders_for_bars[i + 1]) / 2)
        bar_widths.append(borders_for_bars[i + 1] - borders_for_bars[i])

    for border in segmentation_borders:
        ax.axvline(x=border, ymin=0, ymax=0.5, linestyle='--', color='C4')

    all_borders = [1] + segmentation_borders + [ChromPos.chrs[chr] + 1]
    for i in range(len(all_borders) - 1):
        if BADs[i]:
            ax.axhline(y=BADs[i],
                       xmin=all_borders[i] / ChromPos.chrs[chr],
                       xmax=all_borders[i + 1] / ChromPos.chrs[chr],
                       linewidth=BAD_lw, color=BAD_color_1,
                       solid_capstyle='butt')

    ax.scatter(x=chr_snps['pos'], y=list(chr_snps['AD']), c=chr_snps['cov'], cmap='BuPu', s=2, vmin=10, vmax=30)
    ax.set_xlim(0, ChromPos.chrs[chr])
    ax.set_ylim(y_min, y_max)
    ax.grid(which='major', axis='both')
    ax.set_xticklabels([])
    ax.set_yticks(list(range(1, int(y_max) + 1)))
    ax.text(0.99, 0.95, 'K562 {}'.format(chr),
            horizontalalignment='right',
            verticalalignment='top',
            transform=ax.transAxes)
    ax.set_ylabel('AD')
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="10%", pad=0.05)
    cax.get_yaxis().set_ticks([])
    cax.set_xlim(1, ChromPos.chrs[chr])
    cax.set_ylim(0, 1)
    cax.bar(bar_positions, [1] * len(bar_positions), bar_widths, color=reduced_bar_colors, linewidth=0)

cax.set_xlabel('Chromosome position, bp')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)

ax = axs[-1]
ax.plot([0, 0], [0, 0], color=BAD_color_1, label='Estimated BAD')
ax.legend(loc='upper left')

ax = fig.add_axes([0.07, 0.87, 0.2, 0.03])
cmap = 'BuPu'
norm = mpl.colors.Normalize(vmin=10, vmax=30)
cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                               norm=norm,
                               orientation='horizontal')

plt.savefig(os.path.expanduser('~/AC_4/Figure_AS_4_stage2.svg'), dpi=300)
plt.savefig(os.path.expanduser('~/AC_4/Figure_AS_4_stage2.png'), dpi=300)
# plt.show()
plt.close(fig)

#COSMIC comparison
cosmic = pd.read_csv(cosmic_name, low_memory=False)
fig, (*axs,) = plt.subplots(len(chrs), 1)
fig.tight_layout(pad=1.5)
plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))

for chr, ax in zip(chrs, axs):
    chr_ploidy = ploidy[ploidy['#chr'] == chr]
    chr_snps = snps[snps['chr'] == chr].copy()
    chr_snps['AD'] = chr_snps['AD'].apply(lambda y: y_max - delta_y if y > y_max else y)
    chr_cosmic = cosmic.loc[
        (cosmic['#sample_name'] == cnv_line) &
        ('chr' + cosmic['chr'] == chr) &
        (cosmic['minorCN'] != 0)
    ].copy()
    chr_cosmic['chr'] = 'chr' + chr_cosmic['chr']
    chr_cosmic['startpos'] = chr_cosmic['startpos'].astype(int)
    chr_cosmic['endpos'] = chr_cosmic['endpos'].astype(int)

    bar_positions = []
    bar_widths = []
    bar_colors = []
    vd = 1 / 500

    BADs = []

    borders_to_draw = []
    segmentation_borders = []
    last_end = 1
    for index, (pl_chr, start, end, BAD, *values) in chr_ploidy.iterrows():
        if start != last_end:
            if last_end == 1:
                borders_to_draw += [start - ChromPos.chrs[chr] * vd]
                segmentation_borders += [start - ChromPos.chrs[chr] * vd]
            else:
                borders_to_draw += [last_end + ChromPos.chrs[chr] * vd, start - ChromPos.chrs[chr] * vd]
                segmentation_borders += [last_end + ChromPos.chrs[chr] * vd, start - ChromPos.chrs[chr] * vd]
            bar_colors.append('#AAAAAA')
            BADs.append(None)
        else:
            if last_end != 1:
                segmentation_borders += [last_end]
        last_end = end
        bar_colors.append('C2')
        BADs.append(BAD)
    if last_end != ChromPos.chrs[chr] + 1:
        borders_to_draw += [last_end + ChromPos.chrs[chr] * vd]
        segmentation_borders += [last_end + ChromPos.chrs[chr] * vd]
        bar_colors.append('#AAAAAA')
        BADs.append(None)

    reduced_bar_colors = []
    for i, color in enumerate(bar_colors):
        if i == 0 or bar_colors[i - 1] != color:
            reduced_bar_colors.append(color)

    borders_for_bars = [1] + borders_to_draw + [ChromPos.chrs[chr] + 1]
    for i in range(len(borders_for_bars) - 1):
        bar_positions.append((borders_for_bars[i] + borders_for_bars[i + 1]) / 2)
        bar_widths.append(borders_for_bars[i + 1] - borders_for_bars[i])

    for border in segmentation_borders:
        ax.axvline(x=border, ymin=0, ymax=0.5, linestyle='--', color='C4')

    all_borders = [1] + segmentation_borders + [ChromPos.chrs[chr] + 1]
    for i in range(len(all_borders) - 1):
        if BADs[i]:
            ax.axhline(y=BADs[i],
                       xmin=all_borders[i] / ChromPos.chrs[chr],
                       xmax=all_borders[i + 1] / ChromPos.chrs[chr],
                       linewidth=BAD_lw, color=BAD_color,
                       solid_capstyle='butt')

    # cosmic
    cosmic_bar_colors = []
    vd = 1 / 500
    COSMIC_BADs = []
    cosmic_borders = []
    last_end = 1
    for index, (sample_name, sample_id, SNPstart, SNPend, chr,
                startpos, endpos, chr_37, start_37, end_37, minorCN, totalCN) in chr_cosmic.iterrows():
        if startpos - last_end >= ChromPos.chrs[chr] * vd * 2:
            if last_end == 1:
                cosmic_borders += [startpos - ChromPos.chrs[chr] * vd]
            else:
                cosmic_borders += [last_end + ChromPos.chrs[chr] * vd, startpos - ChromPos.chrs[chr] * vd]
            cosmic_bar_colors.append('#AAAAAA')
            COSMIC_BADs.append(None)
        else:
            if last_end != 1:
                cosmic_borders += [last_end]
        last_end = endpos
        cosmic_bar_colors.append('C2')
        COSMIC_BADs.append((totalCN - minorCN) / minorCN)
    if last_end != ChromPos.chrs[chr] + 1:
        cosmic_borders += [last_end + ChromPos.chrs[chr] * vd]
        cosmic_bar_colors.append('#AAAAAA')
        COSMIC_BADs.append(None)

    all_cosmic_borders = [1] + cosmic_borders + [ChromPos.chrs[chr] + 1]
    for i in range(len(all_cosmic_borders) - 1):
        if COSMIC_BADs[i]:
            ax.axhline(y=COSMIC_BADs[i],
                       xmin=all_cosmic_borders[i] / ChromPos.chrs[chr],
                       xmax=all_cosmic_borders[i + 1] / ChromPos.chrs[chr],
                       linewidth=COSMIC_lw, color=COSMIC_color, snap=False, ms=0, mew=0,
                       solid_capstyle='butt')

    ax.scatter(x=chr_snps['pos'], y=list(chr_snps['AD']), c=chr_snps['cov'], cmap='BuPu', s=2, vmin=10, vmax=30)
    ax.set_xlim(0, ChromPos.chrs[chr])
    ax.set_ylim(y_min, y_max)
    ax.grid(which='major', axis='both')
    ax.set_xticklabels([])
    ax.set_yticks(list(range(1, int(y_max) + 1)))
    ax.text(0.99, 0.95, 'K562 {}'.format(chr),
            horizontalalignment='right',
            verticalalignment='top',
            transform=ax.transAxes)
    ax.set_ylabel('AD')
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="10%", pad=0.05)
    cax.get_yaxis().set_ticks([])
    cax.set_xlim(1, ChromPos.chrs[chr])
    cax.set_ylim(0, 1)
    cax.bar(bar_positions, [1] * len(bar_positions), bar_widths, color=reduced_bar_colors, linewidth=0)

cax.set_xlabel('Chromosome position, bp')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)

ax = axs[-1]
ax.plot([0, 0], [0, 0], color=BAD_color, label='Estimated BAD')
ax.plot([0, 0], [0, 0], color=COSMIC_color, label='COSMIC BAD')
ax.legend(loc='upper left')

ax = fig.add_axes([0.07, 0.87, 0.2, 0.03])
cmap = 'BuPu'
norm = mpl.colors.Normalize(vmin=10, vmax=30)
cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                               norm=norm,
                               orientation='horizontal')

plt.savefig(os.path.expanduser('~/AC_4/Figure_AS_4_stage2_cosmic.svg'), dpi=300)
plt.savefig(os.path.expanduser('~/AC_4/Figure_AS_4_stage2_cosmic.png'), dpi=300)
plt.show()
plt.close(fig)