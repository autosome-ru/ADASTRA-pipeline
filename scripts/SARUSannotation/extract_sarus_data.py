import numpy as np
import pandas as pd
from collections import deque
from scripts.HELPERS.helpers import ChromPos
from scripts.HELPERS.paths_for_components import FA, results_path
from scripts.HELPERS.paths import get_tf_sarus_path, get_result_table_path

num = dict(zip(['A', 'a', 'C', 'c', 'G', 'g', 'T', 't'], [1, 1, 2, 2, 3, 3, 4, 4]))
nuc = dict(zip([1, 2, 3, 4, 5], ['a', 'c', 'g', 't', 'N']))


def main(tf_name, motif_length, opened_df=None):
    if motif_length is None:
        exit(0)
    tf_path = get_result_table_path('TF', tf_name)
    out_path = get_tf_sarus_path(tf_name, 'fasta')
    positions = dict()
    gen = dict()
    if opened_df is not None:
        tf_df = opened_df
    else:
        tf_df = pd.read_table(tf_path)
    if tf_df.empty:
        exit(0)
    for chr_name in ChromPos.chromosomes:
        positions[chr_name] = deque()
        gen[chr_name] = np.zeros(ChromPos.chromosomes[chr_name], np.int8)
    for index, row in tf_df.iterrows():
        positions[row['#chr']].append(row['pos'])

    skip_to_next_chr = False
    pos = 0
    l_const = 10000000
    with open(FA, 'r') as fasta:
        for line in fasta:
            if line[0] == '>':
                chromosome = line[1:-1]
                if chromosome not in ChromPos.chromosomes:
                    skip_to_next_chr = True
                    continue
                else:
                    skip_to_next_chr = False
                pos = 0
                count = 0
                next_2 = 0

                next_start = 0
                next_end = 0

                if positions[chromosome]:
                    next_1 = positions[chromosome].popleft()
                    next_start = next_1 - (motif_length - 1)
                    next_end = next_1 + (motif_length - 1)
                    if positions[chromosome]:
                        next_2 = positions[chromosome].popleft()
                        while next_end >= next_2 - (motif_length - 1):
                            next_end = next_2 + (motif_length - 1)
                            if positions[chromosome]:
                                next_2 = positions[chromosome].popleft()
                            else:
                                break
                else:
                    skip_to_next_chr = True
            else:
                if skip_to_next_chr:
                    continue
                if pos < next_start - 110:
                    pos += 100

                    if pos // l_const > count:
                        count = pos // l_const

                    continue
                for sym in line.strip():
                    pos += 1

                    if pos // l_const > count:
                        count = pos // l_const

                    if pos > next_end:
                        if positions[chromosome]:
                            next_1 = next_2
                            next_start = next_1 - (motif_length - 1)
                            next_end = next_1 + (motif_length - 1)
                            if positions[chromosome]:
                                next_2 = positions[chromosome].popleft()
                                while next_end >= next_2 - (motif_length - 1):

                                    next_end = next_2 + (motif_length - 1)
                                    if positions[chromosome]:
                                        next_2 = positions[chromosome].popleft()
                                    else:
                                        break
                        else:
                            if next_2:
                                next_1 = next_2
                                next_2 = 0
                                next_start = next_1 - (motif_length - 1)
                                next_end = next_1 + (motif_length - 1)
                            else:
                                skip_to_next_chr = True
                                break
                    if next_start <= pos <= next_end:
                        gen[chromosome][pos] = num.get(sym.lower(), 5)
    with open(out_path, 'w') as out:
        for index, row in tf_df.iterrows():
            chromosome = row['#chr']
            pos = row['pos']
            R = row['ref']
            A = row['alt']
            ID = row['ID'] + ";" + A
            if gen[chromosome][pos] == 0:
                continue

            assert R.lower() == nuc[gen[chromosome][pos]]
            if pos < 25:
                print('Fuck...')
                exit(1)
            left_tail = ''.join([nuc[gen[chromosome][pos - (motif_length - 1) + i]]
                                 for i in range((motif_length - 1))])
            right_tail = ''.join([nuc[gen[chromosome][pos + 1 + i]]
                                  for i in range((motif_length - 1))])
            out.write('>' + ID + '_ref' + '\n')
            out.write(left_tail + R + right_tail + '\n')

            out.write('>' + ID + '_alt' + '\n')
            out.write(left_tail + A + right_tail + '\n')
