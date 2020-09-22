import json
import sys
import os.path

sys.path.insert(1, "/home/abramov/segmentationValidation/ADASTRA-pipeline")
from scripts.HELPERS.paths_for_components import parallel_parameters_path, ploidy_dict_path
from scripts.HELPERS.helpers import read_synonims

out_path = parallel_parameters_path + 'PE_parameters.cfg'

cosmic_names, cgh_names = read_synonims()

if __name__ == "__main__":
    with open(ploidy_dict_path, 'r') as read_file:
        d = json.loads(read_file.readline())
    keys = sorted(d.keys())
    with open(out_path, 'w') as file:
        for key in keys:
            if not cosmic_names.get(key.split('!')[0]):
                continue
            is_empty = True
            for value in d[key]:
                if os.path.isfile(value):
                    is_empty = False
            if is_empty:
                continue
            file.write(key + '\n')
