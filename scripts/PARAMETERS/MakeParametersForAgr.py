import json
import sys

parameters_path = '/home/abramov/PARAMETERS/'
parallel_parameters_path = '/home/abramov/ParallelParameters/'


if __name__ == "__main__":
    expected_args = {"CL": "TF", "TF": "CL"}
    what_for = sys.argv[1]  # "TF" or "CL" arguments are expected
    if what_for not in expected_args:
        raise ValueError('{} not in CL, TF'.format(what_for))
    with open(parameters_path + what_for + '_DICT.json', 'r') as read_file:
        d = json.loads(read_file.readline())

    keys = sorted(d.keys())
    with open(parallel_parameters_path + 'Agr_parameters.cfg', 'w') as file:
        for key in keys:
            file.write(key + '\n')
