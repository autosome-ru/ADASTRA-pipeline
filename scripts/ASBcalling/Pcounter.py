import json
import sys
from scipy.stats import binom_test
import os.path

ploidy_path = "/home/abramov/PloidyForHotFix/"
parameters_path = "/home/abramov/PARAMETERS/"


class ChromPos:
    def __init__(self, chr, pos):
        self.chr = chr
        self.pos = pos
    
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


def count_p(x, n, p, alternative):
    pv = (binom_test(x, n, p, alternative) + binom_test(x, n, 1-p, alternative))/2
    return pv


def unpack(line):
    line = line.split()
    chr = line[0]
    pos = int(line[1])
    ID = line[2]
    ref = line[3]
    alt = line[4]
    Q = float(line[7])
    (ref_c, alt_c, GQ, in_macs, in_sissrs, in_cpics, in_gem) = map(int, line[5:7]+line[8:13])
    callers = in_macs + in_sissrs + in_cpics + in_gem
    return chr, pos, ID, ref, alt, ref_c, alt_c, Q, GQ, in_macs, in_sissrs, in_cpics, in_gem, callers


def pack(values):
    return '\t'.join(map(str, values)) + '\n'


def create_ploidy(string):
    path = ploidy_path + "Corrected-1,5/" + string + "_ploidy.tsv"
    return path


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

with open(parameters_path + "CELL_LINES.json", "r") as read_file:
    d = json.loads(read_file.readline())
    rev_d = make_reverse_dict(d)

ploidy_file = rev_d.get(key, None)
if ploidy_file is None:
    print("No ploidy found")
    ploidy = None
else:
    ploidy = create_ploidy(ploidy_file)
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
        with open(table_annotated, 'r') as file, open(output, 'w') as out:
            current = 0
            for line in file:
                if line[0] == '#':
                    continue
                
                chr, pos, ID, ref, alt, ref_c, alt_c, Q, GQ, in_macs, in_sissrs, in_cpics, in_gem, callers = unpack(
                                                                                                                line)
                
                # ploidy annotation
                chrom, start, end, ploidy, dip_qual, lq, rq, seg_c = segments[current]
                start = int(start)
                end = int(end)
                cur_st = ChromPos(chrom, start)
                cur_ed = ChromPos(chrom, end)
                now = ChromPos(chr, pos)
                while now >= cur_ed and current+1 < len(segments):
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
                    #p_value counting
                    p = 1 / (float(ploidy) + 1)
                    n = ref_c + alt_c
                    p_ref = count_p(ref_c, n, p, 'greater')
                    p_alt = count_p(ref_c, n, p, 'less')
                else:
                    p_ref = '.'
                    p_alt = '.'
                
                out.write(pack([chr, pos, ID, ref, alt, ref_c, alt_c, Q, GQ, in_macs, in_sissrs, in_cpics,
                                in_gem, callers, ploidy, dip_qual, lq, rq, seg_c, p_ref, p_alt]))
