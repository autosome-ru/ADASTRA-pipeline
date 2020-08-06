import json
import sys
from scripts.HELPERS.paths import create_path_from_master_list
from scripts.HELPERS.paths_for_components import configs_path, master_list_path
from scripts.HELPERS.helpers import check_if_in_expected_args, remove_punctuation


def makedict(what_for):
    d = dict()
    check_if_in_expected_args(what_for)
    with open(master_list_path, "r") as m:
        master = m.readlines()
    for line in master:
        if line[0] == "#":
            continue
        ln = line.split("\t")
        path = create_path_from_master_list(ln, for_what="p-value_table")
        if what_for == "TF":
            try:
                d[ln[1]].append(path)
            except KeyError:
                d[ln[1]] = [path]
        if what_for == "CL":
            cell_line = remove_punctuation(ln[4])
            try:
                d[cell_line].append(path)
            except KeyError:
                d[cell_line] = [path]
    with open(configs_path + what_for + "_DICT.json", "w") as write_file:
        json.dump(d, write_file)
    print("Dictionary Saved")


indicator = sys.argv[1]
makedict(indicator)
