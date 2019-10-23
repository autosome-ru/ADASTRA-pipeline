import json
import sys
from scipy.stats import binom_test
import os.path

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import ploidy_dict_path, create_ploidy_path_function
from scripts.HELPERS.helpers import ChromPos, callers_names, unpack, pack


def count_p(x, n, p, alternative):
    pv = (binom_test(x, n, p, alternative) + binom_test(x, n, 1 - p, alternative)) / 2
    return pv


def make_reverse_dict(dictionary):
    new_dict = {}
    for key in dictionary:
        paths = dictionary[key]
        for path in paths:
            if path.split("/")[-3] != "CTRL":
                new_dict[path] = key
    return new_dict


full_path = sys.argv[1]
key = full_path + ".vcf.gz"

with open(ploidy_dict_path, "r") as read_file:
    d = json.loads(read_file.readline())
    rev_d = make_reverse_dict(d)

ploidy_file = rev_d.get(key, None)
if ploidy_file is None:
    print("No ploidy found")
    ploidy = None
else:
    ploidy = create_ploidy_path_function(ploidy_file)
    if os.path.isfile(ploidy):
        table_annotated = full_path + "_table_annotated.txt"
        output = full_path + "_table_p.txt"
        segments = []
        with open(ploidy, 'r') as file:
            for line in file:
                if line[0] == '#':
                    continue
                line = line.split()
                segments.append(line)

        segments = sorted(segments, key=lambda x: int(x[1]))
        segments = sorted(segments, key=lambda x: x[0])
        if len(segments) == 0:
            print('Ploidy file is empty!')
            exit(1)

        print('Now doing', table_annotated, '\n', 'with ploidy file', ploidy_file)
        with open(output, 'w') as out, open(table_annotated, 'r') as file:
            current = 0
            for line in file:
                if line[0] == '#':
                    continue
                chr, pos, ID, ref, alt, ref_c, alt_c, Q, GQ, in_callers = unpack(line, use_in="Pcounter")

                # ploidy annotation
                chrom, start, end, ploidy, dip_qual, lq, rq, seg_c = segments[current]
                start = int(start)
                end = int(end)
                cur_st = ChromPos(chrom, start)
                cur_ed = ChromPos(chrom, end)
                now = ChromPos(chr, pos)
                while now >= cur_ed and current + 1 < len(segments):
                    current += 1
                    chrom, start, end, ploidy, dip_qual, lq, rq, seg_c = segments[current]
                    start = int(start)
                    end = int(end)
                    cur_st = ChromPos(chrom, start)
                    cur_ed = ChromPos(chrom, end)
                    now = ChromPos(chr, pos)

                if now < cur_st or now >= cur_ed:
                    ploidy = 0
                    p_ref = 0
                    p_alt = 0
                elif float(ploidy) != 0:
                    # p_value counting
                    p = 1 / (float(ploidy) + 1)
                    n = ref_c + alt_c
                    p_ref = count_p(ref_c, n, p, 'greater')
                    p_alt = count_p(ref_c, n, p, 'less')
                else:
                    p_ref = '.'
                    p_alt = '.'

                out.write(pack([chr, pos, ID, ref, alt, ref_c, alt_c, Q, GQ] +
                               [in_callers[name] for name in callers_names] +
                               [ploidy, dip_qual, lq, rq, seg_c, p_ref, p_alt]))
