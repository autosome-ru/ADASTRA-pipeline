import sys
import os
import pandas as pd

class ChromPos:
    chr_l = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
         145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
         101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468,
         156040895, 57227415]
    chrs = dict(zip(['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY'], chr_l))
    genome_length = sum(chr_l)

    def __init__(self, chr, pos):
        if chr not in self.chrs:
            raise ValueError("Not in valid chromosomes {}".format(chr))
        self.chr = chr
        self.pos = int(pos)

    def __lt__(self, other):
        if self.chr == other.chr:
            return self.pos < other.pos
        else:
            return self.chr < other.chr

    def __gt__(self, other):
        if self.chr == other.chr:
            return self.pos > other.pos
        else:
            return self.chr > other.chr

    def __le__(self, other):
        if self.chr == other.chr:
            return self.pos <= other.pos
        else:
            return self.chr <= other.chr

    def __ge__(self, other):
        if self.chr == other.chr:
            return self.pos >= other.pos
        else:
            return self.chr >= other.chr

    def __eq__(self, other):
        return (self.chr, self.pos) == (other.chr, other.pos)

    def __ne__(self, other):
        return (self.chr, self.pos) != (other.chr, other.pos)

    def distance(self, other):
        if self.chr != other.chr:
            return float('inf')
        return abs(self.pos - other.pos)

chrs = ChromPos.chrs.keys()
cosmic_name = os.path.expanduser('~/PARAMETERS/COSMIC_copy_number.sorted.tsv')
ploidy_folder = os.path.expanduser('~/PloidyForRelease/CAIC/')
cnv_line = 'K-562'

cosmic = pd.read_table(cosmic_name, low_memory=False)
cosmic.columns = ['#sample_name', 'chr', 'startpos', 'endpos', 'minorCN', 'totalCN']

cosmics = {}
for chr in chrs:
    chr_cosmic = cosmic.loc[
            (cosmic['#sample_name'] == cnv_line) &
            (cosmic['chr'] == chr) &
            (cosmic['minorCN'] != 0)
        ].copy()
    chr_cosmic['chr'] = chr_cosmic['chr']
    chr_cosmic['startpos'] = chr_cosmic['startpos'].astype(int)
    chr_cosmic['endpos'] = chr_cosmic['endpos'].astype(int)
    cosmics[chr] = chr_cosmic

COSMIC_covers_tr = 0.3
match_tr = 0.95
BAD_segments_tr = 2
different_BADs_tr = 2

with open(os.path.expanduser('~/CHERRIES.txt'), 'w') as out:
    out.write('\t'.join(['file', 'chr', 'COSMIC_covers', 'BAD_match', 'number_of_COSMIC_segments', 'number_of_BAD_segments', 'different_BADs']) + '\n')
    for file in os.listdir(ploidy_folder):
        print(file)
        ploidy = pd.read_table(ploidy_folder + file)
        for chr in chrs:
            chr_ploidy = ploidy[ploidy['#chr'] == chr]

            chromosome_length = ChromPos.chrs[chr]

            total = 0
            intersection = 0

            different_BADs = set()

            number_of_BAD_segments = len(chr_ploidy.index)
            number_of_COSMIC_segments = 0

            BAD_segments_iterator = iter(chr_ploidy.iterrows())
            try:
                index, (pl_chr, start, end, BAD, *values) = next(BAD_segments_iterator)
            except StopIteration:
                continue

            for index, (sample_name, chrom, startpos, endpos, minorCN, totalCN) in cosmics[chr].iterrows():
                total += endpos - startpos
                number_of_COSMIC_segments += 1
                while end < startpos:
                    try:
                        index, (pl_chr, start, end, BAD, *values) = next(BAD_segments_iterator)
                    except StopIteration:
                        break
                while end < endpos:
                    if BAD == totalCN / minorCN - 1:
                        intersection += end - max(start, startpos)
                        different_BADs.add(BAD)
                    try:
                        index, (pl_chr, start, end, BAD, *values) = next(BAD_segments_iterator)
                    except StopIteration:
                        break
                if start < endpos:
                    if BAD == totalCN / minorCN - 1:
                        intersection += endpos - max(start, startpos)
                        different_BADs.add(BAD)

            COSMIC_covers = total / chromosome_length
            if total == 0:
                BAD_match = 0
            else:
                BAD_match = intersection / total

            # print(COSMIC_covers, BAD_match, number_of_COSMIC_segments, number_of_BAD_segments)

            if COSMIC_covers >= COSMIC_covers_tr and BAD_match >= match_tr and number_of_BAD_segments >= BAD_segments_tr and len(different_BADs) >= different_BADs_tr:
                out.write('\t'.join(map(str, [file, chr, COSMIC_covers, BAD_match, number_of_COSMIC_segments, number_of_BAD_segments, different_BADs])) + '\n')
