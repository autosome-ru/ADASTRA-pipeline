import json
import sys
import os.path

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import parameters_path, parallel_parameters_path

out_path = parallel_parameters_path + 'Agr_parameters.cfg'

if __name__ == "__main__":
    expected_args = {"CL": "TF", "TF": "CL"}
    what_for = sys.argv[1]  # "TF" or "CL" arguments are expected
    if what_for not in expected_args:
        raise ValueError('{} not in CL, TF'.format(what_for))
    with open(parameters_path + what_for + '_DICT.json', 'r') as read_file:
        d = json.loads(read_file.readline())

    keys = sorted(d.keys())
    with open(out_path, 'w') as file:
        for key in keys:
            is_empty = True
            for value in d[key]:
                if os.path.isfile(value):
                    is_empty = False
            if is_empty:
                continue
            file.write(key + '\n')
