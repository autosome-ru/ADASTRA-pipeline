import json
import sys

if __name__ == "__main__":
    TF_json_path = '/home/abramov/ASBcalling/TF_DICT.json'
    CL_json_path = '/home/abramov/ASBcalling/CL_DICT.json'
    out_path = '/home/abramov/ParallelParameters/PvC_parameters.cfg'

    what_for = sys.argv[1]
    if what_for == "TF":
        with open(TF_json_path, 'r') as read_file:
            d = json.loads(read_file.readline())
    elif what_for == "CL":
        with open(CL_json_path, 'r') as read_file:
            d = json.loads(read_file.readline())

    keys = sorted(d.keys())
    with open(out_path, 'w') as file:
        for key in keys:
            file.write(key + '\n')
