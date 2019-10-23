import json
import sys

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import ploidy_dict_path, parallel_parameters_path

out_path = parallel_parameters_path + 'PE_parameters.cfg'

if __name__ == "__main__":
    with open(ploidy_dict_path, 'r') as read_file:
        d = json.loads(read_file.readline())

    keys = sorted(d.keys())
    with open(out_path, 'w') as file:
        for key in keys:
            file.write(key + '\n')
