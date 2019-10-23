from statistics import mean, median_grouped
from scipy.stats import kendalltau

callers_names = ['macs', 'sissrs', 'cpics', 'gem']

chr_l = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
         145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
         101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468,
         156040895, 57227415]

Nucleotides = {'A', 'T', 'G', 'C'}


class ChromPos:
    chrs = dict(zip(['chr' + str(i) for i in range(1, 22)] + ['chrX', 'chrY'], chr_l))
    
    def __init__(self, chr, pos):
        assert chr in self.chrs
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


def unpack_segments(line):
    if line[0] == '#':
        return [None] * 3
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
        self.snp_args = None
        self.seg_args = None
        self.snp_coordinate = None
        self.segment_start = None
        self.segment_end = None
        self.has_segments = True
        
    def __iter__(self):
        return self
    
    def return_snp(self, intersect):
        return (self.snp_coordinate.chr, self.snp_coordinate.pos, self.snp_args) \
               + (int(intersect),) * self.write_intersect \
               + tuple(arg * intersect for arg in self.seg_args) * self.write_segment_args
    
    def get_next_snp(self):
        try:
            snp_chr, pos, *self.snp_args = self.unpack_snp_function(next(self.snps))
            if snp_chr is None:
                self.get_next_snp()
                return
            self.snp_coordinate = ChromPos(snp_chr, pos)
        except ValueError:
            raise StopIteration
    
    def get_next_segment(self):
        try:
            seg_chr, start_pos, end_pos, *self.seg_args = self.unpack_segments_function(next(self.segments))
            if seg_chr is None:
                self.get_next_segment()
                return
            self.segment_start = ChromPos(seg_chr, start_pos)
            self.segment_end = ChromPos(seg_chr, end_pos)
        except (StopIteration, ValueError):
            self.has_segments = False
        
    def __next__(self):
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
            self.get_next_snp()
            return x


def make_dict_from_vcf(vcf, vcf_dict):
    for line in vcf:
        if line[0] == '#':
            continue
        line = line.split()
        chr = line[0]
        if chr not in ChromPos.chrs:
            raise ValueError('{} not in valid chrs'.format(chr))
        pos = int(line[1])
        if not len(line[3]) == 1 or not len(line[4]) == 1:
            continue
        if line[3] not in Nucleotides or line[4] not in Nucleotides:
            continue
        Inf = line[-1].split(':')
        R = int(Inf[1].split(',')[0])
        if Inf[1].split(",")[1] == "":
            print(line)
            print(vcf)
        A = int(Inf[1].split(',')[1])
        if min(R, A) < 3:
            continue
        GT = Inf[0]
        if GT != '0/1':
            continue
        ID = line[2]
        REF = line[3]
        ALT = line[4]
        try:
            prev_value = vcf_dict[(chr, pos, ID, REF, ALT)]
            vcf_dict[(chr, pos, ID, REF, ALT)] = (R + prev_value[0], A + prev_value[1])
        except KeyError:
            vcf_dict[(chr, pos, ID, REF, ALT)] = (R, A)


def unpack(line, use_in):
    if use_in == "Pcounter" and line[0] == '#':
        return [None] * 3
    line = line.strip().split('\t')
    chr = line[0]
    pos = int(line[1])
    ID = line[2]
    ref = line[3]
    alt = line[4]
    ref_c, alt_c = map(int, line[5:7])
    if use_in == "PloidyEstimation":
        return chr, pos, ID, ref, alt, ref_c, alt_c
    Q = float(line[7])
    GQ = int(line[8])
    difference = len(callers_names)

    peaks = map(int, line[9:9 + difference])
    in_callers = dict(zip(callers_names, [peaks]))
    if use_in == "Pcounter":
        return chr, pos, ID, ref, alt, ref_c, alt_c, Q, GQ, in_callers

    if use_in == "Aggregation":
        dip_qual, lq, rq, seg_c = map(int, line[10 + difference:14 + difference])
        ploidy = float(line[9 + difference])
        if line[19] == '.':
            p_ref = '.'
            p_alt = '.'
        else:
            p_ref, p_alt = map(float, line[14 + difference:16 + difference])
        return chr, pos, ID, ref, alt, ref_c, alt_c, Q, GQ, in_callers, ploidy, \
               dip_qual, lq, rq, seg_c, p_ref, p_alt

    raise ValueError('{} not in Aggregation, P-value, PloidyEstimation options for function usage'.format(use_in))


def pack(values):
    return '\t'.join(map(str, values)) + '\n'
