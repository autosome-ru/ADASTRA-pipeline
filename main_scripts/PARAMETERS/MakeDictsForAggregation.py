import os.path
import json
import sys

alignments_path = "/home/abramov/Alignments/"
master_path = "/home/abramov/PLOIDYcalling/"
dict_path = "/home/abramov/ASBcalling/"


def createpath(line):
    return alignments_path + "EXP/" + line[1] + "/" + line[0] + "/" + line[3] + "_table_p.txt"


def makedicts():
    d = dict()
    with open(master_path + "Master-lines.tsv", "r") as m:
        master = m.readlines()
    for line in master:
        ln = line.split()
        path = createpath(ln)
        if os.path.isfile(path):
            if what_for == "TF":
                try:
                    d[ln[1]].append(path)
                except KeyError:
                    d[ln[1]] = [path]
            elif what_for == "CL":
                #FIXME
                try:
                    d[ln[1]].append(path)
                except KeyError:
                    d[ln[1]] = [path]
    with open(dict_path + what_for + "_DICT.json", "w") as write_file:
        json.dump(d, write_file)
    print("Dictionaries Saved")


makedicts()
