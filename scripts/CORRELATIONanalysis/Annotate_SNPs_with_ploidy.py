import os
import sys
import json

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import ploidy_dict_path, ploidy_path
from scripts.HELPERS.helpers import Intersection

scripts_path = '/home/abramov/ASB-Project/scripts'
Correlation_path = '/home/abramov/Correlation/'


def unpack_ploidy_segments(line):
    line = line.strip().split('\t')
    return [line[0], int(line[1]), int(line[2]), float(line[3]), int(line[4]), int(line[7])]


def unpack_snps(line):
    line = line.strip().split('\t')
    return [line[0], int(line[1]), int(line[5]), int(line[6])]


names = []
with open(scripts_path + '/CORRELATIONanalysis/synonims.tsv', 'r') as syn:
    for line in syn:
        line = line.strip('\n').split('\t')
        if line[1]:
            GTRD_name = line[0].replace('(', '').replace(')', '').replace(' ', '_')
            COSMIC_name = line[1]
            names.append(GTRD_name)

count = dict()
for name in names:
    count[name] = 0

with open(ploidy_dict_path, 'r') as file:
    cl = json.loads(file.readline().strip())

modes = []
for file_name in sorted(os.listdir(ploidy_path)):
    if os.path.isdir(ploidy_path + file_name):
        modes.append(file_name)
        
file_name = sys.argv[1]

assert os.path.isfile(ploidy_path + file_name)

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
    table_path = ploidy_path + file_name
    for mode in modes:
        if not os.path.isdir(Correlation_path + mode + '_tables/'):
            os.mkdir(Correlation_path + mode + '_tables')
        ploidy_file_path = ploidy_path + mode + '/' + name + '!' + lab + '_ploidy.tsv'
        out_path = Correlation_path + mode + '_tables/' + name + '_' + lab.replace('_', '-') + '.tsv'
        print(out_path)

        with open(table_path, 'r') as table, open(ploidy_file_path, 'r') as ploidy, open(out_path, 'w') as out:
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
                segments.append(Segment(line[0], int(line[1]), int(line[2]), float(line[3]), int(line[4]), int(line[7])))
                segments = sorted(segments, key=lambda x: x.start)

            out.write('##' + str(len(result)) + '!' + str(datasetsn) + '!' + str(sum(1 for segment in segments if segment.value >= 1)) + '!'+ lab+ '!' +'>'.join(al_list))
            out.write('\t'.join(['#chr', 'pos', 'ref', 'alt', 'ploidy', 'qual', 'segn']) + '\n')
            for line in result:
                out.write('\t'.join(map(str, line)) + '\n')

