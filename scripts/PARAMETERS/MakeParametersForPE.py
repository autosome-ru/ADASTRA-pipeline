import json
import sys
import os.path

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import ploidy_dict_path, parallel_parameters_path
from scripts.HELPERS.helpers import read_synonims

out_path = parallel_parameters_path + 'PE_parameters.cfg'

cell_lines = ['K-562', 'MCF7', 'HCT-116']
test_names = [
    "K562_myelogenous_leukemia!_labs_michael-snyder___biosamples_ENCBS357NWO__ploidy.tsv",
    "K562_myelogenous_leukemia!_labs_michael-snyder___biosamples_ENCBS603CUX__ploidy.tsv",
    "K562_myelogenous_leukemia!_labs_michael-snyder___biosamples_ENCBS906KIP__ploidy.tsv",
    "K562_myelogenous_leukemia!_labs_michael-snyder___biosamples_ENCBS389PVA__ploidy.tsv",
    "K562_myelogenous_leukemia!labs_michael-snyder___biosamples_ENCBS392AAA__ploidy.tsv",
    "K562_myelogenous_leukemia!_labs_bradley-bernstein___biosamples_ENCBS639AAA__ploidy.tsv",
    "K562_myelogenous_leukemia!_labs_richard-myers___biosamples_ENCBS577JGM__ploidy.tsv",
              ]

if __name__ == "__main__":
    with open(ploidy_dict_path, 'r') as read_file:
        d = json.loads(read_file.readline())
    cosmic_names, _ = read_synonims()
    keys = sorted(d.keys())
    with open(out_path, 'w') as file:
        for key in keys:
            print(key)
            if key not in test_names:
                continue
            name = key.split('!')[0]
            if cosmic_names.get(name, '') not in cell_lines:
                continue
            is_empty = True
            for value in d[key]:
                if os.path.isfile(value):
                    is_empty = False
            if is_empty:
                continue
            file.write(key + '\n')
