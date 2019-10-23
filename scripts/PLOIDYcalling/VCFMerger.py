import gzip
import os
import sys
import json

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import ploidy_path, ploidy_dict_path
from scripts.HELPERS.helpers import make_dict_from_vcf


def merge_vcfs(out_file_name, in_files):
    vcf_dict = dict()
    for file in in_files:
        with gzip.open(file, 'rt') as vcf:
            make_dict_from_vcf(vcf, vcf_dict)
    
    vcf_keys = list(vcf_dict.keys())
    vcf_keys.sort(key=lambda cords: cords[1])
    vcf_keys.sort(key=lambda cords: cords[0])
    
    with open(out_file_name, 'w') as out:
        for (chr, pos, ID, REF, ALT) in vcf_keys:
            (R, A) = vcf_dict[(chr, pos, ID, REF, ALT)]
            out.write('\t'.join(map(str, [chr, pos, ID, REF, ALT, R, A])) + '\n')


if __name__ == '__main__':
    with open(ploidy_dict_path, 'r') as read_file:
        d = json.loads(read_file.readline())
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
