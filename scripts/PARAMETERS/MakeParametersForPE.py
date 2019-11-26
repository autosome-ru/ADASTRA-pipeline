import json
import sys
import os.path

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import ploidy_dict_path, parallel_parameters_path
from scripts.HELPERS.helpers import read_synonims

out_path = parallel_parameters_path + 'PE_parameters.cfg'

cell_lines = ['K-562', 'MCF7', 'HCT-116']
test_names = [
    # "K562_myelogenous_leukemia!_labs_michael-snyder___biosamples_ENCBS357NWO_", it's too big
    "K562_myelogenous_leukemia!_labs_michael-snyder___biosamples_ENCBS603CUX_",
    "K562_myelogenous_leukemia!_labs_michael-snyder___biosamples_ENCBS906KIP_",
    "K562_myelogenous_leukemia!_labs_michael-snyder___biosamples_ENCBS389PVA_",
    "K562_myelogenous_leukemia!_labs_michael-snyder___biosamples_ENCBS392AAA_",
    "K562_myelogenous_leukemia!_labs_bradley-bernstein___biosamples_ENCBS639AAA_",
    "K562_myelogenous_leukemia!_labs_richard-myers___biosamples_ENCBS577JGM_",
              ]

if __name__ == "__main__":
    with open(ploidy_dict_path, 'r') as read_file:
        d = json.loads(read_file.readline())
    #cosmic_names, _ = read_synonims()
    keys = sorted(d.keys())
    with open(out_path, 'w') as file:
        for key in keys:
            # name = key.split('!')[0]
            # if cosmic_names.get(name, '') not in cell_lines:
            #     continue
            is_empty = True
            for value in d[key]:
                if os.path.isfile(value):
                    is_empty = False
            if is_empty:
                continue
            file.write(key + '\n')
