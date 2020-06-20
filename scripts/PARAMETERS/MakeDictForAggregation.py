import json
import sys

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import create_path_from_GTRD_function
from scripts.HELPERS.paths_for_components import parameters_path, GTRD_slice_path
from scripts.HELPERS.helpers import check_if_in_expected_args, remove_punctuation


def makedict(what_for):
    d = dict()
    check_if_in_expected_args(what_for)
    with open(GTRD_slice_path, "r") as m:
        master = m.readlines()
    for line in master:
        if line[0] == "#":
            continue
        ln = line.split("\t")
        path = create_path_from_GTRD_function(ln, for_what="p-value_table")
        if what_for == "TF":
            try:
                d[ln[1]].append(path)
            except KeyError:
                d[ln[1]] = [path]
        if what_for == "CL":
            cell_line = ln[4].replace("(", "").replace(")", "").replace(" ", "_").replace("/", "_")
            try:
                d[cell_line].append(path)
            except KeyError:
                d[cell_line] = [path]
    with open(parameters_path + what_for + "_DICT.json", "w") as write_file:
        json.dump(d, write_file)
    print("Dictionary Saved")


indicator = sys.argv[1]
makedict(indicator)
