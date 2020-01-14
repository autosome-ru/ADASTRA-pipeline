import sys
import os.path
from statistics import median_grouped
from scipy import stats
import numpy as np
import json
import statsmodels.stats.multitest
import pandas as pd
import string
from shutil import copyfile


sys.path.insert(1, "/home/abramov/ASB-Project")

from scripts.HELPERS.paths_for_components import parameters_path, results_path,\
    GTRD_slice_path
from scripts.HELPERS.helpers import callers_names, unpack, pack, check_if_in_expected_args


def makedict():
    d = dict()
    check_if_in_expected_args(what_for)
    with open(GTRD_slice_path, "r") as m:
        master = m.readlines()
    for line in master:
        if line[0] == "#":
            continue
        ln = line.split("\t")
        cell_line = remove_punctuation(ln[4]).replace(" ", "_")
        d[cell_line] = ln[4]
    return d


def remove_punctuation(x):
    table = str.maketrans({key: "_" for key in string.punctuation})
    return x.translate(table)


what_for = 'CL'
convert_cl = makedict()
with open(parameters_path + 'CONVERT_CL_NAMES.json') as out:
    json.dump(convert_cl, out)