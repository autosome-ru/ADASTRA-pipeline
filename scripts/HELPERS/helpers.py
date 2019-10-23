from statistics import mean, median_grouped
from scipy.stats import kendalltau

callers_names = ['macs', 'sissrs', 'cpics', 'gem']

class ChromPos:
    chrs = {'chrX', 'chrY'}
    for i in range(1, 22):
        chrs.add('chr' + str(i))
    
    def __init__(self, chr, pos):
        assert chr in self.chrs
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
    
    def distance(self, other):
        if self.chr != other.chr:
            return float('inf')
        return abs(self.pos - other.pos)


Nucleotides = {'A', 'T', 'G', 'C'}


def make_dict_from_vcf(vcf, vcf_dict):
    strange = 0
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
        NAME = line[2]
        REF = line[3]
        ALT = line[4]
        prev_value = vcf_dict.get((chr, pos), None)
        if prev_value:
            if NAME != prev_value[2] or REF != prev_value[3] or ALT != prev_value[4]:
                strange += 1
                continue
            vcf_dict[(chr, pos)] = (R + prev_value[0], A + prev_value[1], NAME, REF, ALT)
        else:
            vcf_dict[(chr, pos)] = (R, A, NAME, REF, ALT)
    
    return strange


class GObject:
    def __init__(self, chr, pos, value, qual, snpn):
        self.chr_pos = ChromPos(chr, pos)
        self.value = value
        self.qual = qual
        self.snpn = snpn


class Segment:
    def __init__(self, chr, st, ed, value):
        assert (ed > st)
        self.start = ChromPos(chr, st)
        self.end = ChromPos(chr, ed)
        self.value = value
    
    def length(self):
        return self.end.pos - self.start.pos


class ObjectTable:
    def __init__(self, objects=None):
        if objects:
            self.objects = objects
            self.sort_items()
        else:
            self.objects = []
    
    def sort_items(self):
        self.objects = sorted(self.objects, key=lambda x: x.chr_pos)
    
    def add_object(self, obj):
        self.objects.append(obj)
    
    def add_new_objects(self, objects):
        self.objects += objects
        self.sort_items()
    
    def merge_with_object_table(self, object_table):
        pass
    
    def merge_with_segment_table(self, segment_table):
        # print('merging SNP')
        if len(segment_table.segments) == 0:
            return []
        other_current_idx = 0
        other_max_idx = len(segment_table.segments) - 1
        other_current_start = segment_table.segments[other_current_idx].start
        other_current_end = segment_table.segments[other_current_idx].end
        result = []  # Object list
        for obj in self.objects:
            while obj.chr_pos >= other_current_end and other_current_idx + 1 <= other_max_idx:
                other_current_idx += 1
                other_current_start = segment_table.segments[other_current_idx].start
                other_current_end = segment_table.segments[other_current_idx].end
            if obj.chr_pos < other_current_start or obj.chr_pos >= other_current_end:
                continue
            result.append((
                obj.chr_pos.chr,
                obj.chr_pos.pos,
                obj.value,
                obj.qual,
                obj.snpn,
                segment_table.segments[other_current_idx].value,
                segment_table.segments[other_current_idx].start.distance(
                    segment_table.segments[other_current_idx].end),
            ))
        # result = sorted(result, key=lambda x: ChromPos(x[0], x[1]))
        # print(len(self.objects), len(result))
        return result
    
    def print_merged_to_file(self, segment_table, filename):
        with open(filename, 'w') as out:
            out.write('\t'.join(['#chr', 'pos', 'value', 'qual', 'snpn', 'segment_value', 'segment_length']) + "\n")
            result = self.merge_with_segment_table(segment_table)
            for line in result:
                out.write('\t'.join(map(str, line)) + '\n')
            return result
    
    @staticmethod
    def correlation_of_merged(result):
        # print('calc cor of SNP')
        idx = ['chr', 'pos', 'value', 'qual', 'snpn', 'segment_value', 'segment_length'].index('segment_value')
        values1 = []
        values2 = []
        for line in result:
            values1.append(line[2])
            values2.append(line[idx])
        return kendalltau(values1, values2)[0]


class SegmentTable:
    def __init__(self, segments=None):
        if segments:
            self.segments = segments
            self.sort_items()
        else:
            self.segments = []
    
    def sort_items(self):
        self.segments = sorted(self.segments, key=lambda x: x.start)
    
    def add_segment(self, segment):
        self.segments.append(segment)
    
    def add_new_segments(self, segments):
        self.segments += segments
        self.sort_items()
    
    def merge_with_object_table(self, object_table):
        # print('merging COSMIC')
        if len(object_table.objects) == 0:
            return []
        other_current_idx = 0
        other_max_idx = len(object_table.objects) - 1
        other_current_chrpos = object_table.objects[other_current_idx].chr_pos
        result = []  # Segment list
        for segment in self.segments:
            vals = []
            while other_current_chrpos < segment.start and other_current_idx + 1 <= other_max_idx:
                other_current_idx += 1
                other_current_chrpos = object_table.objects[other_current_idx].chr_pos
            while other_current_chrpos <= segment.end and other_current_idx + 1 <= other_max_idx:
                vals.append(object_table.objects[other_current_idx].value)
                other_current_idx += 1
                other_current_chrpos = object_table.objects[other_current_idx].chr_pos
            if not vals:
                continue
            assert len(vals) > 0
            # print(segment.value)
            result.append((
                segment.start.chr,
                segment.start.pos,
                segment.end.pos,
                segment.value,
                len(vals),
                mean(vals),
                median_grouped(vals),
            ))
        # result = sorted(result, key=lambda x: ChromPos(x[0], x[1]))
        # print(len(self.segments), len(result))
        return result
    
    def merge_with_segment_table(self, segment_table):
        pass
    
    def print_merged_to_file(self, object_table, filename):
        with open(filename, 'w') as out:
            out.write('\t'.join(['#chr', 'start', 'end', 'value', 'mean', 'med', 'count']))
            for line in self.merge_with_object_table(object_table):
                out.write('\t'.join(map(str, line)) + '\n')
    
    @staticmethod
    def correlation_of_merged(method, result):
        # print('calc cor of COSMIC')
        idx = ['chr', 'start', 'end', 'value', 'count', 'mean', 'med'].index(method)
        values1 = []
        values2 = []
        for line in result:
            values1.append(line[3])
            values2.append(line[idx])
        return kendalltau(values1, values2)[0]


class Reader:
    CGH_path = ''
    SNP_path = ''
    Cosmic_path = ''
    synonims_path = ''
    
    def read_Cosmic(self, name, mode='normal'):
        with open(self.Cosmic_path, 'r') as file:
            result = SegmentTable()
            for line in file:
                if line[0] == '#':
                    continue
                line = line.strip().split(',')
                # if int(line[4]) in {4,6,8} or line[3] == '0': continue
                if line[0] != name:
                    continue
                if 'chr' + line[4] not in ChromPos.chrs:
                    continue
                if int(line[10]) == 0:
                    continue
                
                if mode == 'normal':
                    value = int(line[11]) / int(line[10]) - 1
                elif mode == 'total':
                    value = int(line[11])
                else:
                    raise ValueError(mode)
                
                result.add_segment(Segment('chr' + line[4], int(line[5]), int(line[6]), value))
            if not result.segments:
                raise KeyError(name)
            # result.sort_items()
            return result
    
    def read_SNPs(self, method='normal'):
        with open(self.SNP_path, 'r') as file:
            result = ObjectTable()
            for line in file:
                if line[0] == '#':
                    idx = line[2:line.rfind('#')].split('!')
                    aligns = idx[4]
                    lab = idx[3]
                    segsegs = idx[2]
                    datas = idx[1]
                    idx = idx[0]
                    if aligns:
                        aligns = ','.join(aligns.split('>'))
                    else:
                        aligns = ''
                    continue
                line = line.split()
                if line[0] not in ChromPos.chrs:
                    continue
                if method == 'normal':
                    if line[4] == 0:
                        continue
                    result.add_object(GObject(line[0], int(line[1]), float(line[4]), int(line[5]), int(line[6])))
                elif method == 'naive':
                    ref = int(line[2])
                    alt = int(line[3])
                    if min(ref, alt) == 0:
                        continue
                    result.add_object(GObject(line[0], int(line[1]), max(ref, alt) / min(ref, alt) - 1, 10000, 10000))
                else:
                    raise KeyError(method)
            # result.sort_items()
            return idx, datas, lab, result, aligns, segsegs
    
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
        idx = cgnames.index(cgh_name) + 3
        N = 0
        with open(self.CGH_path, 'r') as file:
            result = ObjectTable()
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
                N += 1
                result.add_object(GObject(chr, pos, value, 100, 100))
            # result.sort_items()
            return N, result
    
    def read_synonims(self):
        cosmic_names = dict()
        cgh_names = dict()
        with open(self.synonims_path, 'r') as file:
            for line in file:
                line = line.strip('\n').split('\t')
                if line[1] and line[2]:
                    name = line[0].replace(')', '').replace('(', '').replace(' ', '_')
                    cosmic_names[name] = line[1]
                    cgh_names[name] = line[2]
        return cosmic_names, cgh_names
