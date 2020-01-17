import sys
import numpy as np
from collections import deque

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.helpers import ChromPos

fasta = open(sys.argv[2], 'r')
table = open(sys.argv[1], 'r')

num = dict(zip(['A', 'a', 'C', 'c', 'G', 'g', 'T', 't'], [1, 1, 2, 2, 3, 3, 4, 4]))
nuc = dict(zip([1, 2, 3, 4, 5], ['a', 'c', 'g', 't', 'N']))

positions = dict()
gen = dict()

for chr_name in ChromPos.chrs:
    positions[chr_name] = deque()
    gen[chr_name] = np.zeros(ChromPos.chrs[chr_name], np.int8)

for line in table:
    if line[0] == '#':
        continue
    line = line.split()
    positions[line[0]].append(int(line[1]))
table.seek(0)

ok = False
for key in positions:
    if positions[key] != deque():
        ok = True

if not ok:
    print('No snps!')
    exit(0)

out = open(sys.argv[3], 'w')
motive_length = int(sys.argv[4])

p = 0
l_const = 10000000
skip_to_next_chr = False
for line in fasta:
    if line[0] == '>':
        chr = line[1:-1]
        print(chr)
        if chr not in ChromPos.chrs:
            skip_to_next_chr = True
            continue
        else:
            skip_to_next_chr = False
        p = 0

        count = 0

        next_1 = 0
        next_2 = 0

        if positions[chr]:
            next_1 = positions[chr].popleft()
            next_start = next_1 - (motive_length - 1)
            next_end = next_1 + (motive_length - 1)
            # print(next_start, next_end)
            if positions[chr]:
                next_2 = positions[chr].popleft()
                # print(next_2)
                while next_end >= next_2 - (motive_length - 1):
                    next_end = next_2 + (motive_length - 1)
                    if positions[chr]:
                        next_2 = positions[chr].popleft()
                    else:
                        break
        else:
            skip_to_next_chr = True
        print(next_start, next_end)
    else:
        if skip_to_next_chr:
            continue
        if p < next_start - 110:
            p += 100

            if p // l_const > count:
                count = p // l_const
            # print(count*l, '+100', next_start, next_end)

            continue
        for sym in line.strip():
            p += 1

            if p // l_const > count:
                count = p // l_const
            # print(count*l, '+1', next_start, next_end)

            if p > next_end:
                if positions[chr]:
                    next_1 = next_2
                    next_start = next_1 - (motive_length - 1)
                    next_end = next_1 + (motive_length - 1)
                    if positions[chr]:
                        next_2 = positions[chr].popleft()
                        while next_end >= next_2 - (motive_length - 1):
                            # print('a-haa!', p, next_start, next_end, next_2-(motive_length - 1))
                            next_end = next_2 + (motive_length - 1)
                            if positions[chr]:
                                next_2 = positions[chr].popleft()
                            else:
                                break
                else:
                    if next_2:
                        next_1 = next_2
                        next_2 = 0
                        next_start = next_1 - (motive_length - 1)
                        next_end = next_1 + (motive_length - 1)
                    else:
                        skip_to_next_chr = True
                        break
            if next_start <= p <= next_end:
                gen[chr][p] = num.get(sym.lower(), 5)

for line in table:
    assert type(line) == str
    if line[0] == '#':
        continue
    line = line.split()
    chr = line[0]
    try:
        p = int(line[1])
    except ValueError:
        continue
    R = line[3]
    A = line[4]
    # print(chr, p, gen[chr][p])
    if gen[chr][p] == 0:
        continue
    # print(chr, p, gen[chr][p])
    if R.lower() != nuc[gen[chr][p]]:
        print('Vse ploho', chr, line[1])
    out.write(chr + '_' + line[1] + ' ' + ''.join(
        [nuc[gen[chr][p - (motive_length - 1) + i]] for i in range((motive_length - 1))]) +
              '[' + R + '/' + A + ']' + ''.join([nuc[gen[chr][p + 1 + i]] for i in range((motive_length - 1))]) + '\n')
