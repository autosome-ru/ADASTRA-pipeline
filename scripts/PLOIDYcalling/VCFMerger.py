import gzip
import os
import sys
import json

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import parameters_path, ploidy_path

Nucleotides = {'A', 'T', 'G', 'C'}


def make_dict_from_vcf(vcf, vcf_dict):
    strange = 0
    for line in vcf:
        if line[0] == '#':
            continue
        line = line.split()
        if len(line) != 10:
            print('Shit happens')
            print('\t'.join(line))
            continue
        chr = line[0]
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


def merge_vcfs(out_file_name, in_files):
    vcf_dict = dict()
    strange = 0
    for file in in_files:
        with gzip.open(file, 'rt') as vcf:
            strange += make_dict_from_vcf(vcf, vcf_dict)
    
    vcf_keys = list(vcf_dict.keys())
    vcf_keys.sort(key=lambda cords: cords[1])
    vcf_keys.sort(key=lambda cords: cords[0])
    
    with open(out_file_name, 'w') as out:
        for (chr, pos) in vcf_keys:
            (R, A, NAME, REF, ALT) = vcf_dict.get((chr, pos))
            out.write('\t'.join(map(str, [chr, pos, NAME, REF, ALT, R, A])) + '\n')
    
    print('Skipped {0} mismatched SNPS'.format(strange))


if __name__ == '__main__':
    with open(parameters_path + 'CELL_LINES.json', 'r') as read_file:
        d = json.loads(read_file.readline())
    keys = sorted(d.keys())
    
    key = sys.argv[1]
    print(key)
    
    arr = []
    
    for path in d[key]:
        if os.path.isfile(path):
            arr.append(path)
        else:
            continue
    if not arr:
        sys.exit(0)
    
    out_file = ploidy_path + key + ".tsv"
    print(arr)
    
    merge_vcfs(out_file, arr)
