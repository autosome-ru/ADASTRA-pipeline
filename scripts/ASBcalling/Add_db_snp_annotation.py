import json
import sys
import numpy as np
import pandas as pd
import gzip


sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import create_ploidy_path_function
from scripts.HELPERS.paths_for_components import parameters_path, ploidy_dict_path
from scripts.HELPERS.helpers import callers_names, unpack, pack, Intersection

if __name__ == '__main__':
    with open(parameters_path + 'RS_list.tsv') as f:
        rs_array = set(json.loads(f.readline()))
    list_of_tags = ['NSF', 'NSM', 'NSN', 'REF', 'SYN', 'U3', 'U5', 'ASS', 'DSS', 'INT', 'R3', 'R5']
    header = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    with gzip.open('/home/abramov/REFERENCE/00-common_all.vcf.gz', 'rt') as rf, \
            open(parameters_path + 'rs_info.tsv', 'w') as out:
        out.write(pack(['ID'] + list_of_tags))
        counter = 0
        for line in rf:
            if counter % 20000 == 0:
                print(counter)
            counter += 1
            if line[0] == '#':
                continue
            line = line.strip('\n').split('\t')
            ID = line[2]
            if ID not in rs_array:
                continue
            info = line[-1].split(';')
            out.write(pack([ID] + [1 if tag in info else 0 for tag in list_of_tags]))



