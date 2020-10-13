import string
import numpy as np
import os
from .paths_for_components import master_list_path, configs_path, synonyms_path

from .paths import create_neg_bin_weights_path_function

callers_names = ['macs', 'sissrs', 'cpics', 'gem', 'macs2', 'macs2-nomodel']

chr_l = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
         145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
         101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468,
         156040895, 57227415]

Nucleotides = {'A', 'T', 'G', 'C'}
expected_args = {"CL": "TF", "TF": "CL"}

states = [1, 3/2, 2, 3, 4, 5]

master_list_header = '#EXP	TF_UNIPROT_ID	ANTIBODY	TREATMENT	SPECIE	CELL_ID	CELLS	EXP_TYPE	CONTROL	READS	ALIGNS	PEAKS	GEO	ENCODE	WG_ENCODE	READS_ALIGNED   DOWNLOAD_PATH'

dtype_dict = {name: str if name != 'READS_ALIGNED' else np.float_ for name in master_list_header.split('\t')}


class ChromPos:
    chromosomes = dict(zip(['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY'], chr_l))
    genome_length = sum(chr_l)

    def __init__(self, chromosome, pos):
        if chromosome not in self.chromosomes:
            raise ValueError("Not in valid chromosomes {}".format(chromosome))
        self.chromosome = chromosome
        self.pos = int(pos)

    def __lt__(self, other):
        if self.chromosome == other.chromosome:
            return self.pos < other.pos
        else:
            return self.chromosome < other.chromosome

    def __gt__(self, other):
        if self.chromosome == other.chromosome:
            return self.pos > other.pos
        else:
            return self.chromosome > other.chromosome

    def __le__(self, other):
        if self.chromosome == other.chromosome:
            return self.pos <= other.pos
        else:
            return self.chromosome <= other.chromosome

    def __ge__(self, other):
        if self.chromosome == other.chromosome:
            return self.pos >= other.pos
        else:
            return self.chromosome >= other.chromosome

    def __eq__(self, other):
        return (self.chromosome, self.pos) == (other.chromosome, other.pos)

    def __ne__(self, other):
        return (self.chromosome, self.pos) != (other.chromosome, other.pos)

    def distance(self, other):
        if self.chromosome != other.chromosome:
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
        return [self.snp_coordinate.chromosome, self.snp_coordinate.pos] + self.snp_args \
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
        if chr not in ChromPos.chromosomes:
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


def make_list_from_vcf(vcf, file_name=None, filter_no_rs=False):
    vcf_list = []
    for line in vcf:
        if line[0] == '#':
            continue
        line = line.split()
        chr = line[0]
        if chr not in ChromPos.chromosomes:
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
        if file_name is not None:
            vcf_list.append((chr, pos, ID, REF, ALT, R, A, file_name))
        else:
            vcf_list.append((chr, pos, ID, REF, ALT, R, A, file_name))
    return vcf_list


def make_list_from_vcf_without_filter(vcf):
    vcf_list = []
    for line in vcf:
        if line[0] == '#':
            continue
        line = line.split()
        chr = line[0]
        if chr not in ChromPos.chromosomes:
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


def check_if_in_expected_args(what_for):
    if what_for not in expected_args:
        raise ValueError('{} not in CL, TF'.format(what_for))


def remove_punctuation(x):
    table = str.maketrans({key: "_" for key in string.punctuation if key not in {'-', '+'}})
    return x.translate(table).replace(" ", "_")


def read_weights():
    r = {}
    w = {}
    gof = {}
    for fixed_allele in ('ref', 'alt'):
        r[fixed_allele] = {}
        w[fixed_allele] = {}
        gof[fixed_allele] = {}
        for BAD in states:
            precalc_params_path = create_neg_bin_weights_path_function(fixed_allele, BAD)
            coefs_array = np.load(precalc_params_path)
            r[fixed_allele][BAD] = coefs_array[:, 0]
            w[fixed_allele][BAD] = coefs_array[:, 1]
            gof[fixed_allele][BAD] = coefs_array[:, 3]
            good_gof_list = [x for x in range(len(gof[fixed_allele][BAD])) if gof[fixed_allele][BAD][x] > 0.05]
            if len(good_gof_list) > 0:
                first_bad_gof = min(x for x in range(len(gof[fixed_allele][BAD])) if gof[fixed_allele][BAD][x] > 0.05)
            else:
                first_bad_gof = 0
            gof[fixed_allele][BAD][first_bad_gof:] = 1
            r[fixed_allele][BAD][first_bad_gof:] = 0
            w[fixed_allele][BAD][first_bad_gof:] = 1
    return r, w, gof


class UnpackBadSegments:
    counter = None

    def __init__(self, counter):
        UnpackBadSegments.counter = counter

    def unpackBADSegments(self, line, states):
        if line[0] == '#':
            return [''] * (len(line.strip().split('\t')) - len(states) + 1 + (1 if UnpackBadSegments.counter is not None else 0))
        line = line.strip().split('\t')

        if UnpackBadSegments.counter is not None:
            UnpackBadSegments.counter += 1
            return [line[0], int(line[1]), int(line[2]), float(line[3]), UnpackBadSegments.counter] + \
                   [dict(zip(states, line[4: 4 + len(states)]))] + line[(4 + len(states)):]
        else:
            return [line[0], int(line[1]), int(line[2]), float(line[3])] + \
                   [dict(zip(states, line[4: 4 + len(states)]))] + line[(4 + len(states)):]


def get_states(states_sign):
    if states_sign == 'all_but_1.33_2.5':
        states = [1, 2, 3, 4, 5, 6, 1.5]
    elif states_sign == '123456':
        states = [1, 2, 3, 4, 5, 6]
    elif states_sign == '12345':
        states = [1, 2, 3, 4, 5]
    elif states_sign == 'all_but_1.33':
        states = [1, 2, 3, 4, 5, 1.5, 6, 2.5]
    elif states_sign == 'all_but_2.5':
        states = [1, 2, 3, 4, 5, 1.5, 6, 4/3]
    elif states_sign == 'all':
        states = [1, 2, 3, 4, 5, 1.5, 6, 4/3, 2.5]
    elif states_sign == 'all_5':
        states = [1, 2, 3, 4, 5, 1.5]
    else:
        states = [1, 2, 3, 4, 5, 1.5]
    return sorted(states)


class CorrelationReader:
    CGH_path = ''
    SNP_path = ''
    states = []

    def read_SNPs(self, method='normal'):
        with open(self.SNP_path, 'r') as file:
            result = []
            uniq_segments_count = 0
            sum_cov = 0
            previous_segment = []
            for line in file:
                line = line.strip()
                if line[0] == '#':
                    split_header = line[1:].strip().split('@')
                    datasets_number = split_header[0]
                    lab = split_header[1]
                    aligns = split_header[2]
                    if aligns:
                        aligns = ','.join(aligns.split(','))
                    else:
                        aligns = ''
                    continue
                line = line.split("\t")
                if line[0] not in ChromPos.chromosomes:
                    continue
                current_segment = [float(line[4]),
                                   dict(zip(self.states, list(map(float, line[5:5 + len(self.states)]))))]
                if previous_segment != current_segment:
                    uniq_segments_count += 1
                    sum_cov += int(line[6 + len(self.states)])
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
                                   dict(zip(self.states, [10000] * len(self.states)))])
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
                if chr not in ChromPos.chromosomes:
                    continue
                pos = (int(line[1]) + int(line[2])) // 2
                try:
                    value = 2 ** (1 + float(line[idx]))
                except ValueError:
                    continue
                result.append([chr, pos, value, dict(zip(states, [10000] * len(states)))])
            # result.sort_items()
            return result


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


if __name__ == "__main__":
    import pandas as pd
    import json
    with open(os.path.join(configs_path, 'CONVERT_CL_NAMES.json'), 'w') as o:
        d_to_write = {}
        d = pd.read_table(master_list_path, dtype=dtype_dict)
        for key in d["CELLS"].tolist():
            d_to_write[key] = remove_punctuation(key)
        json.dump(d_to_write, o)
