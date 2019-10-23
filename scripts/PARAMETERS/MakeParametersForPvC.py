import json
import os.path
import sys

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import parallel_parameters_path, ploidy_dict_path

out_path = parallel_parameters_path + 'PvC_parameters.cfg'


def make_reverse_dict(dictionary):
    new_dict = {}
    for key in dictionary:
        paths = dictionary[key]
        for path in paths:
            if path.split("/")[-3] != "CTRL":
                new_dict[path] = key
    return new_dict


if __name__ == "__main__":
    with open(ploidy_dict_path, 'r') as read_file:
        d = json.loads(read_file.readline())
        rev_d = make_reverse_dict(d)
    keys = sorted(rev_d.keys())

    with open(out_path, 'w') as file:
        for k in keys:
            full_path = ".".join(k.split(".")[:-2])
            if os.path.isfile(k) and os.path.isfile(full_path + "_table_annotated.txt"):
                file.write(full_path + '\n')
