import json
import sys

parameters_path = '/home/abramov/PARAMETERS/'
parallel_parameters_path = '/home/abramov/ParallelParameters/'


if __name__ == "__main__":
    what_for = sys.argv[1]
    if what_for == "TF":
        with open(parameters_path + 'TF_DICT.json', 'r') as read_file:
            d = json.loads(read_file.readline())
    elif what_for == "CL":
        with open(parameters_path + 'CL_DICT.json', 'r') as read_file:
            d = json.loads(read_file.readline())

    keys = sorted(d.keys())
    with open(parallel_parameters_path + 'Agr_parameters.cfg', 'w') as file:
        for key in keys:
            file.write(key + '\n')
