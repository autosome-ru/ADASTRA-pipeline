import sys
import os
import pandas as pd

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.helpers import ChromPos

chrs = ChromPos.chrs.keys()
cosmic_name = os.path.expanduser('~/Documents/ASB/Cell_lines/cell_lines_copy_number.csv')
cnv_line = 'K-562'

cosmic = pd.read_csv(cosmic_name, low_memory=False)

cosmics = {}
for chr in chrs:
    chr_cosmic = cosmic.loc[
            (cosmic['#sample_name'] == cnv_line) &
            ('chr' + cosmic['chr'] == chr) &
            (cosmic['minorCN'] != 0)
        ].copy()
    chr_cosmic['chr'] = 'chr' + chr_cosmic['chr']
    chr_cosmic['startpos'] = chr_cosmic['startpos'].astype(int)
    chr_cosmic['endpos'] = chr_cosmic['endpos'].astype(int)
    cosmics[chr] = chr_cosmic

COSMIC_covers_tr = 0.3
match_tr = 0.95
BAD_segments_tr = 2
different_BADs_tr = 2

with open(os.path.expanduser('~/CHERRIES.txt'), 'w') as out:
    out.write('\t'.join(['file', 'chr', 'COSMIC_covers', 'BAD_match', 'number_of_COSMIC_segments', 'number_of_BAD_segments', 'different_BADs']) + '\n')
    for file in os.listdir(os.path.expanduser('~/K562_BAD_Segments')):
        print(file)
        ploidy = pd.read_table(os.path.expanduser('~/K562_BAD_Segments/') + file)
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

            for index, (sample_name, sample_id, SNPstart, SNPend, _,
                        startpos, endpos, chr_37, start_37, end_37, minorCN, totalCN) in cosmics[chr].iterrows():
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
