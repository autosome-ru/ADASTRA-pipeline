import json
import sys

alignments_path = "/home/abramov/Alignments/"
master_path = "/home/abramov/PLOIDYcalling/"
dict_path = "/home/abramov/ASBcalling/"


def createpath(line):
    return alignments_path + "EXP/" + line[1] + "/" + line[0] + "/" + line[3] + "_table_p.txt"


def makedict(what_for):
    d = dict()
    with open(master_path + "Master-lines.tsv", "r") as m:
        master = m.readlines()
    for line in master:
        if line[0] == "#":
            continue
        ln = line.split("\t")
        path = createpath(ln)
        if what_for == "TF":
            if d[ln[1]] == "None":
                continue
            try:
                d[ln[1]].append(path)
            except KeyError:
                d[ln[1]] = [path]
        elif what_for == "CL":
            cell_line = ln[4].replace("(", "").replace(")", "").replace(" ", "_").replace("/", "_")
            try:
                d[cell_line].append(path)
            except KeyError:
                d[cell_line] = [path]
    with open(dict_path + what_for + "_DICT.json", "w") as write_file:
        json.dump(d, write_file)
    print("Dictionary Saved")


indicator = sys.argv[1]
makedict(indicator)
