import os
import sys
from statistics import mean, median_grouped
import json


class ChromPos:
    chrs = {'chrX', 'chrY'}
    for i in range(1, 23):
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


class GObject:
    def __init__(self, chr, pos, ref, alt):
        self.chr_pos = ChromPos(chr, pos)
        self.ref = ref
        self.alt = alt


class Segment:
    def __init__(self, chr, st, ed, value, qual, segn):
        assert (ed > st)
        self.start = ChromPos(chr, st)
        self.end = ChromPos(chr, ed)
        self.value = value
        self.qual = qual
        self.segn = segn
    
    def length(self):
        return self.end.pos - self.start.pos


syn_path = '/home/abramov/ASB-Project/main_scripts/CORRELATIONanalysis/synonims.tsv'
JSON_path = '/home/abramov/PLOIDYcalling/CELL_LINES.json'
Ploidy_path = '/home/abramov/Ploidy/'
Correlation_path = '/home/abramov/Correlation'
names = []
with open(syn_path, 'r') as syn:
    for line in syn:
        line = line.strip('\n').split('\t')
        if line[1] and line[2]:
            GTRD_name = line[0].replace('(', '').replace(')', '').replace(' ', '_')
            COSMIC_name = line[1]
            names.append(GTRD_name)

count = dict()
for name in names:
    count[name] = 0

print(names)

with open(JSON_path, 'r') as file:
    cl = json.loads(file.readline().strip())

modes = []
for file_name in sorted(os.listdir(Ploidy_path)):
    if os.path.isdir(Ploidy_path + file_name):
        modes.append(file_name)
        
file_name = sys.argv[1]

assert os.path.isfile(Ploidy_path + file_name)

name = file_name.split('!')[0]
lab = file_name.split('!')[1][:-4]
try:
    aligns = list(set(cl[file_name[:-4]]))
    datasetsn = len(aligns)
    al_list = [align[29:-7] for align in aligns]
except KeyError:
    datasetsn = 'nan'
    al_list = []
    print(file_name)
if name in names:
    count[name] = count[name] + 1
    table_path = Ploidy_path + file_name
    for mode in modes:
        if not os.path.isdir(Ploidy_path + mode + '_tables/'):
            os.mkdir(Ploidy_path + mode + '_tables')
        ploidy_path = Ploidy_path + mode + '/' + name + '!' + lab + '_ploidy.tsv'
        out_path = Correlation_path + mode + '_tables/' + name + '_' + lab.replace('_', '-') + '.tsv'
        print(out_path)
    
        with open(table_path, 'r') as table, open(ploidy_path, 'r') as ploidy:
            objects = []
            segments = []
            for line in table:
                if line[0] == '#':
                    continue
                line = line.split()
                objects.append(GObject(line[0], int(line[1]), int(line[5]), int(line[6])))
            for line in ploidy:
                if line[0] == '#':
                    continue
                line = line.split()
                segments.append(Segment(line[0], int(line[1]), int(line[2]), int(line[3]), int(line[4]), int(line[7])))
        
        with open(out_path, 'w') as out:
            other_current_idx = 0
            other_max_idx = len(segments) - 1
            other_current_start = segments[other_current_idx].start
            other_current_end = segments[other_current_idx].end
            result = []  # Object list
            for object in objects:
                while object.chr_pos >= other_current_end and other_current_idx + 1 <= other_max_idx:
                    other_current_idx += 1
                    other_current_start = segments[other_current_idx].start
                    other_current_end = segments[other_current_idx].end
                if object.chr_pos < other_current_start or object.chr_pos >= other_current_end:
                    continue
                result.append((
                    object.chr_pos.chr,
                    object.chr_pos.pos,
                    object.ref,
                    object.alt,
                    segments[other_current_idx].value,
                    segments[other_current_idx].qual,
                    segments[other_current_idx].segn,
                ))
            result = sorted(result, key=lambda x: ChromPos(x[0], x[1]))
            out.write('##' + str(len(result)) + '!' + str(datasetsn) + '!' + str(len(segments)) + '!'+ lab+ '!' +'>'.join(al_list))
            out.write('\t'.join(['#chr', 'pos', 'ref', 'alt', 'ploidy', 'qual', 'segn']) + '\n')
            for line in result:
                out.write('\t'.join(map(str, line)) + '\n')

