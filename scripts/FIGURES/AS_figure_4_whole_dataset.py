import sys
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt, ticker
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
from scipy.stats import kendalltau

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.helpers import ChromPos

sns.set(font_scale=1.2, style="ticks", font="lato", palette=('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2',
                                                             '#D55E00', '#CC79A7'))
# sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})
plt.rcParams['font.weight'] = "medium"
plt.rcParams['axes.labelweight'] = 'medium'
plt.rcParams['figure.titleweight'] = 'medium'
plt.rcParams['axes.titleweight'] = 'medium'
plt.rcParams['figure.figsize'] = 14, 2
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

snps_name = os.path.expanduser('~/Desktop//SRR8652105.table.tsv')
ploidy_name = os.path.expanduser('~/Desktop//SRR8652105.badmap.tsv')


# cosmic_name = os.path.expanduser('~/Desktop/CCLE_copy_number_hg38.sorted.txt')
cosmic_name = os.path.expanduser('~/Desktop/COSMIC_copy_number_MCF7.txt')
cnv_line = 'MCF7_BREAST'

snps = pd.read_table(snps_name, header=None)
snps.columns = ['chr', 'pos', 'ID', 'ref', 'alt', 'ref_c', 'alt_c', 'dataset']
snps = snps[snps['ID'] != '.']
snps['AD'] = snps[['ref_c', 'alt_c']].max(axis=1) / snps[['ref_c', 'alt_c']].min(axis=1)
snps['cov'] = snps['ref_c'] + snps['alt_c']
snps['log_cov'] = np.log10(snps['cov'])

ploidy = pd.read_table(ploidy_name)

# chrs = ('chr10', 'chr17')
chrs = ChromPos.chrs.keys()
BAD_color = '#0072B266'
BAD_color_1 = '#0072B2CC'
COSMIC_color = '#D55E00'
BAD_lw = 8
COSMIC_lw = 3.5
y_min = 0.8
y_max = 6
delta_y = 0.05

def unpack_segments(line):
    if isinstance(line, (list, tuple)):
        return line
    if line[0] == '#':
        return [''] * len(line.strip().split('\t'))
    return line.strip().split('\t')

class Intersection:
    def __init__(self, snps, segments, write_segment_args=False, write_intersect=False,
                 unpack_segments_function=unpack_segments, unpack_snp_function=lambda x: x):
        self.snps = iter(snps)
        self.segments = iter(segments)
        self.unpack_snp_function = unpack_snp_function
        self.unpack_segments_function = unpack_segments_function
        self.write_segment_args = write_segment_args
        self.write_intersect = write_intersect
        self.snp_args = []
        self.seg_args = []
        self.snp_coordinate = None
        self.segment_start = None
        self.segment_end = None
        self.has_segments = True
        self.has_snps = True

    def __iter__(self):
        return self

    def return_snp(self, intersect):
        return [self.snp_coordinate.chr, self.snp_coordinate.pos] + self.snp_args \
               + [int(intersect)] * self.write_intersect \
               + [((arg if intersect else {}) if isinstance(arg, dict) else arg * intersect)
                  for arg in self.seg_args] * self.write_segment_args

    def get_next_snp(self):
        try:
            snp_chr, pos, *self.snp_args = self.unpack_snp_function(next(self.snps))
            self.snp_coordinate = ChromPos(snp_chr, pos)
        except ValueError:
            self.get_next_snp()

    def get_next_segment(self):
        try:
            seg_chr, start_pos, end_pos, *self.seg_args = self.unpack_segments_function(next(self.segments))
            self.segment_start = ChromPos(seg_chr, start_pos)
            self.segment_end = ChromPos(seg_chr, end_pos)
        except StopIteration:
            self.has_segments = False
        except ValueError:
            self.get_next_segment()

    def __next__(self):
        if not self.has_snps:
            raise StopIteration

        if self.snp_coordinate is None:
            self.get_next_snp()
        if self.segment_start is None:
            self.get_next_segment()

        while self.has_segments and self.snp_coordinate >= self.segment_end:
            self.get_next_segment()

        if self.has_segments and self.snp_coordinate >= self.segment_start:
            x = self.return_snp(True)
            self.get_next_snp()
            return x
        else:
            x = self.return_snp(False)
            try:
                self.get_next_snp()
            except StopIteration:
                self.has_snps = False
            return x

def correlation_with_cosmic():
    cosmic_segments = []
    with open(cosmic_name, 'r') as cosmic_file:
        for line in cosmic_file:
            if line.startswith('#'):
                continue
            sample_name, chr, startpos, endpos, minorCN, totalCN = line.strip('\n').split('\t')
            if sample_name != 'MCF7_BREAST':
                continue
            if int(minorCN) == 0:
                continue
            cosmic_segments.append((chr, int(startpos), int(endpos), (int(totalCN) - int(minorCN))/int(minorCN) ))
    ploidy_segments = []
    with open(ploidy_name, 'r') as ploidy_file:
        for line in ploidy_file:
            if line.startswith('#'):
                continue
            pl_chr, start, end, BAD, *values = line.strip('\n').split('\t')
            if BAD == 0:
                continue
            ploidy_segments.append((pl_chr, int(start), int(end), float(BAD)))
    snps = []
    with open(snps_name, 'r') as snps_file:
        for line in snps_file:
            if line.startswith('#'):
                continue
            chr, pos, ID, ref, alt, ref_c, alt_c, dataset = line.strip('\n').split('\t')
            snps.append((chr, int(pos)))


    pl_snps = []
    for chr, pos, in_intersect, BAD \
            in Intersection(snps, ploidy_segments, write_intersect=True,
                            write_segment_args=True):
        if not in_intersect:
            continue
        pl_snps.append((chr, pos, BAD))


    snp_ploidy = []
    cosm_ploidy = []
    for chr, pos, ploidy, in_intersect, cosmic_ploidy \
            in Intersection(pl_snps, cosmic_segments, write_intersect=True,
                            write_segment_args=True):
        if not in_intersect:
            continue
        snp_ploidy.append(ploidy)
        cosm_ploidy.append(cosmic_ploidy)

    if len(snp_ploidy) != 0:
        kt = kendalltau(snp_ploidy, cosm_ploidy)[0]
        if kt == 'nan':
            return 'NaN'
        return kt
    return 'NaN'

print(correlation_with_cosmic())


#COSMIC comparison
for chr in chrs:
    print(chr)
    cosmic = pd.read_table(cosmic_name, low_memory=False)
    fig, ax = plt.subplots(1, 1)
    fig.tight_layout(pad=1.5)
    plt.gca().xaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))

    chr_ploidy = ploidy[ploidy['#chr'] == chr]
    chr_snps = snps[snps['chr'] == chr].copy()
    chr_snps['AD'] = chr_snps['AD'].apply(lambda y: y_max - delta_y if y > y_max else y)
    chr_cosmic = cosmic.loc[
        (cosmic['#sample_name'] == cnv_line) &
        (cosmic['chr'] == chr) &
        (cosmic['minorCN'] != 0)
    ].copy()
    chr_cosmic['chr'] = chr_cosmic['chr']
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

    if BADs != [None]:
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
    for index, (sample_name, chr, startpos, endpos, minorCN, totalCN) in chr_cosmic.iterrows():
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
    if COSMIC_BADs != [None]:
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
    ax.text(0.99, 0.95, 'MCF7 {}'.format(chr),
            horizontalalignment='right',
            verticalalignment='top',
            transform=ax.transAxes)
    ax.set_ylabel('AD')
    # plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="10%", pad=0.05)
    cax.get_yaxis().set_ticks([])
    cax.set_xlim(1, ChromPos.chrs[chr])
    cax.set_ylim(0, 1)
    cax.bar(bar_positions, [1] * len(bar_positions), bar_widths, color=reduced_bar_colors, linewidth=0)

    if chr == 'chrY':
        cax.set_xlabel('Chromosome position, bp')
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)

    # ax.plot([0, 0], [0, 0], color=BAD_color, label='Estimated BAD')
    # ax.plot([0, 0], [0, 0], color=COSMIC_color, label='COSMIC BAD')
    # ax.legend(loc='upper left')
    #
    # ax = fig.add_axes([0.07, 0.87, 0.2, 0.03])
    # cmap = 'BuPu'
    # norm = mpl.colors.Normalize(vmin=10, vmax=30)
    # cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
    #                                norm=norm,
    #                                orientation='horizontal')

    # plt.savefig(os.path.expanduser('~/AC_4/Figure_AS_4_stage2_cosmic_{}.svg'.format(chr)), dpi=300)
    plt.savefig(os.path.expanduser('~/Desktop/SRR/Figure_AS_4_stage2_cosmic_{}.png'.format(chr)), dpi=300)
    # plt.show()
    plt.close(fig)
