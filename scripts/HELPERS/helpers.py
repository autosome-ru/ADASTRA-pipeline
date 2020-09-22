import sys
import string
import numpy as np

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import make_black_list, create_path_from_GTRD_function
from scripts.HELPERS.paths_for_components import GTRD_slice_path, synonims_path, parameters_path, cl_dict_path

callers_names = ['macs', 'sissrs', 'cpics', 'gem']

chr_l = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
         145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
         101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468,
         156040895, 57227415]

Nucleotides = {'A', 'T', 'G', 'C'}
expected_args = {"CL": "TF", "TF": "CL"}

states = [1, 4/3, 3/2, 2, 5/2, 3, 4, 5, 6]


class ChromPos:
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


def make_dict_from_vcf(vcf, vcf_dict):
    for line in vcf:
        if line[0] == '#':
            continue
        line = line.split()
        chr = line[0]
        if chr not in ChromPos.chrs:
            continue
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
        if min(R, A) < 5:
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


def make_list_from_vcf(vcf, filter_no_rs=False):
    vcf_list = []
    for line in vcf:
        if line[0] == '#':
            continue
        line = line.split()
        chr = line[0]
        if chr not in ChromPos.chrs:
            continue
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
        if min(R, A) < 5:
            continue
        GT = Inf[0]
        if GT != '0/1':
            continue
        ID = line[2]
        if not ID.startswith('rs') and filter_no_rs:
            continue
        REF = line[3]
        ALT = line[4]
        vcf_list.append((chr, pos, ID, REF, ALT, R, A))
    return vcf_list


def make_list_from_vcf_without_filter(vcf):
    vcf_list = []
    for line in vcf:
        if line[0] == '#':
            continue
        line = line.split()
        chr = line[0]
        if chr not in ChromPos.chrs:
            continue
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
        GT = Inf[0]
        if GT != '0/1':
            continue
        ID = line[2]
        REF = line[3]
        ALT = line[4]
        vcf_list.append((chr, pos, ID, REF, ALT, R, A))
    return vcf_list


def unpack(line, use_in):
    if line[0] == '#':
        return []
    line_split = line.strip('\n').split('\t')
    chr = line_split[0]
    pos = int(line_split[1])
    ID = line_split[2]
    ref = line_split[3]
    alt = line_split[4]
    ref_c, alt_c = map(int, line_split[5:7])
    if use_in == "PloidyEstimation":
        return chr, pos, ID, ref, alt, ref_c, alt_c
    repeat = line_split[7]
    difference = len(callers_names)
    peaks = list(map(int, line_split[8:8 + difference]))
    in_callers = dict(zip(callers_names, peaks))
    if use_in == "Pcounter":
        return chr, pos, ID, ref, alt, ref_c, alt_c, repeat, in_callers
    ploidy = float(line_split[8 + difference])
    quals = list(map(float, line_split[9 + difference:9 + difference + len(states)]))
    quals_dict = dict(zip(states, quals))
    difference += len(states)
    seg_c, sum_cov = map(int, line_split[9 + difference:11 + difference])
    p_ref, p_alt = map(float, line_split[11 + difference:13 + difference])
    es_ref = line_split[13 + difference]
    es_alt = line_split[14 + difference]
    if es_ref != '' and es_ref is not None:
        es_ref = float(es_ref)
    else:
        es_ref = None
    if es_alt != '' and es_alt is not None:
        es_alt = float(es_alt)
    else:
        es_alt = None
    if use_in == "Aggregation":
        return chr, pos, ID, ref, alt, ref_c, alt_c, repeat, in_callers, ploidy, quals_dict,\
               seg_c, sum_cov, p_ref, p_alt, es_ref, es_alt

    raise ValueError('{} not in Aggregation, Pcounter, PloidyEstimation options for function usage'.format(use_in))


def pack(values):
    return '\t'.join(map(str, values)) + '\n'


def make_list_for_VCFs(out_path=None, condition_function=lambda x: True):
    # condition function takes path and return boolean
    black_list = make_black_list()
    counted_controls = set()
    condition_to_file_dict = {}
    with open(GTRD_slice_path, "r") as master_list:
        if out_path is not None:
            out = open(out_path, "w")
        for line in master_list:
            if line[0] == "#":
                continue
            split_line = line.strip().split("\t")
            if split_line[0] not in black_list:
                vcf_path = create_path_from_GTRD_function(split_line, for_what="vcf")
                if out_path is not None:
                    if condition_function(vcf_path):
                        out.write(create_line_for_snp_calling(split_line))
                else:
                    condition_to_file_dict[vcf_path] = condition_function(vcf_path)
            if len(split_line) > 10 and split_line[10] not in black_list:
                vcf_path = create_path_from_GTRD_function(split_line, for_what="vcf", ctrl=True)
                if vcf_path in counted_controls:
                    continue
                counted_controls.add(vcf_path)
                if out_path is not None:
                    if condition_function(vcf_path):
                        out.write(create_line_for_snp_calling(split_line, is_ctrl=True))
                else:
                    condition_to_file_dict[vcf_path] = condition_function(vcf_path)
        if out_path is not None:
            out.close()
        else:
            return condition_to_file_dict


def check_if_in_expected_args(what_for):
    if what_for not in expected_args:
        raise ValueError('{} not in CL, TF'.format(what_for))


def read_synonims():
    cosmic_names = dict()
    cgh_names = dict()
    with open(synonims_path, 'r') as file:
        for line in file:
            line = line.strip('\n').split('\t')
            name = remove_punctuation(line[0])
            cosmic_names[name] = line[1]
            cgh_names[name] = line[2]
    return cosmic_names, cgh_names


def remove_punctuation(x):
    table = str.maketrans({key: "_" for key in string.punctuation if key not in {'-', '+'}})
    return x.translate(table).replace(" ", "_")


class CorrelationReader:
    CGH_path = ''
    SNP_path = ''

    def read_SNPs(self, method='normal'):
        with open(self.SNP_path, 'r') as file:
            result = []
            uniq_segments_count = 0
            sum_cov = 0
            previous_segment = []
            for line in file:
                line = line.strip()
                if line[0] == '#':
                    split_header = line[1:].split('!')
                    datasets_number = split_header[0]
                    lab = split_header[1]
                    aligns = split_header[2]
                    if aligns:
                        aligns = ','.join(aligns.split('>'))
                    else:
                        aligns = ''
                    continue
                line = line.split("\t")
                if line[0] not in ChromPos.chrs:
                    continue
                current_segment = [float(line[4]), dict(zip(states, list(map(float, line[5:5 + len(states)]))))]
                if previous_segment != current_segment:
                    uniq_segments_count += 1
                    sum_cov += int(line[6 + len(states)])
                    previous_segment = current_segment
                if method == 'normal':
                    if line[4] == 0:
                        continue
                    result.append([line[0], int(line[1])] + current_segment)
                elif method == 'cover':
                    if line[4] == 0:
                        continue
                    result.append([line[0], int(line[1]), int(line[2]) + int(line[3])] + current_segment)
                elif method == 'naive':
                    ref = int(line[2])
                    alt = int(line[3])
                    if min(ref, alt) == 0:
                        continue
                    result.append([line[0], int(line[1]), max(ref, alt) / min(ref, alt) - 1,
                                   dict(zip(states, [10000] * len(states)))])
                else:
                    raise KeyError(method)

            return datasets_number, lab, result, aligns, uniq_segments_count, sum_cov
    
    def read_CGH(self, cgh_name):
        cgnames = ['BR:MCF7', 'BR:MDA-MB-231', 'BR:HS 578T', 'BR:BT-549', 'BR:T-47D', 'CNS:SF-268', 'CNS:SF-295',
                   'CNS:SF-539', 'CNS:SNB-19', 'CNS:SNB-75', 'CNS:U251', 'CO:COLO 205', 'CO:HCC-2998', 'CO:HCT-116',
                   'CO:HCT-15', 'CO:HT29', 'CO:KM12', 'CO:SW-620', 'LE:CCRF-CEM', 'LE:HL-60(TB)', 'LE:K-562',
                   'LE:MOLT-4', 'LE:RPMI-8226', 'LE:SR', 'ME:LOX IMVI', 'ME:MALME-3M', 'ME:M14', 'ME:SK-MEL-2',
                   'ME:SK-MEL-28', 'ME:SK-MEL-5', 'ME:UACC-257', 'ME:UACC-62', 'ME:MDA-MB-435', 'ME:MDA-N',
                   'LC:A549/ATCC', 'LC:EKVX', 'LC:HOP-62', 'LC:HOP-92', 'LC:NCI-H226', 'LC:NCI-H23', 'LC:NCI-H322M',
                   'LC:NCI-H460', 'LC:NCI-H522', 'OV:IGROV1', 'OV:OVCAR-3', 'OV:OVCAR-4', 'OV:OVCAR-5', 'OV:OVCAR-8',
                   'OV:SK-OV-3', 'OV:NCI/ADR-RES', 'PR:PC-3', 'PR:DU-145', 'RE:786-0', 'RE:A498', 'RE:ACHN',
                   'RE:CAKI-1', 'RE:RXF 393', 'RE:SN12C', 'RE:TK-10', 'RE:UO-31']
        if cgh_name not in cgnames:
            return []
        idx = cgnames.index(cgh_name) + 3
        with open(self.CGH_path, 'r') as file:
            result = []
            for line in file:
                line = line.strip().split('\t')
                chr = line[0]
                if chr not in ChromPos.chrs:
                    continue
                pos = (int(line[1]) + int(line[2])) // 2
                try:
                    value = 2 ** (1 + float(line[idx]))
                except ValueError:
                    continue
                result.append([chr, pos, value, dict(zip(states, [10000]* len(states)))])
            # result.sort_items()
            return result


def correct_cosmic_file(cosmic_file, out_file):
    with open(cosmic_file, 'r') as file, open(out_file, 'w') as out:
        lines = dict()
        for line in file:
            if line[0] == '#':
                continue
            line = line.strip().split(',')
            if not line[0]:
                continue
            if 'chr' + line[4] not in ChromPos.chrs:
                continue
            try:
                lines[line[0]].append(['chr' + line[4], line[5], line[6], line[10], line[11]])
            except KeyError:
                lines[line[0]] = [['chr' + line[4], line[5], line[6], line[10], line[11]]]

        for cell_line, segments in lines.items():
            segments = sorted(segments, key=lambda x: int(x[1]))
            segments = sorted(segments, key=lambda x: x[0])
            for segment in segments:
                out.write(pack([cell_line] + segment))


def create_line_for_snp_calling(split_line, is_ctrl=False):
    if is_ctrl:
        result = [split_line[10]] + ["None", "Homo sapiens"] + split_line[11:15]
        return pack(result)
    else:
        return pack(split_line[:7])


def read_weights():
    r = {}
    w = {}
    gof = {}
    for fixed_allele in ('ref', 'alt'):
        r[fixed_allele] = {}
        w[fixed_allele] = {}
        gof[fixed_allele] = {}
        for BAD in states:
            precalc_params_path = parameters_path + 'NBweights_{}_BAD={:.1f}.npy'.format(fixed_allele, BAD)
            coefs_array = np.load(precalc_params_path)
            r[fixed_allele][BAD] = coefs_array[:, 0]
            w[fixed_allele][BAD] = coefs_array[:, 1]
            gof[fixed_allele][BAD] = coefs_array[:, 3]
            first_bad_gof = min(x for x in range(len(gof[fixed_allele][BAD])) if gof[fixed_allele][BAD][x] > 0.05)
            gof[fixed_allele][BAD][first_bad_gof:] = 1
            r[fixed_allele][BAD][first_bad_gof:] = 0
            w[fixed_allele][BAD][first_bad_gof:] = 1
    return r, w, gof


def unpackBADSegments(line, states):
    if line[0] == '#':
        return [''] * (len(line.strip().split('\t')) - len(states) + 1)
    line = line.strip().split('\t')

    return [line[0], int(line[1]), int(line[2]), float(line[3])] + \
           [dict(zip(states, line[4: 4 + len(states)]))] + line[(4 + len(states)):]


if __name__ == "__main__":
    import pandas as pd
    import json
    with open(parameters_path + 'CONVERT_CL_NAMES.json', 'w') as o:
        d_to_write = {}
        d = pd.read_table(GTRD_slice_path)
        for key in d["cell_title"].tolist():
            d_to_write[key] = remove_punctuation(key)
        json.dump(d_to_write, o)

