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
from scripts.HELPERS.helpers import callers_names, unpack, pack, check_if_in_expected_args, expected_args, states, \
    encode_GTRD_cell_line_name


def makedict():
    d = dict()
    check_if_in_expected_args(what_for)
    with open(GTRD_slice_path, "r") as m:
        master = m.readlines()
    for line in master:
        if line[0] == "#":
            continue
        ln = line.split("\t")
        cell_line = encode_GTRD_cell_line_name(ln[4])
        d[cell_line] = remove_punctuation(ln[4]).replace(" ", "_")
    return d


def remove_punctuation(x):
    table = str.maketrans({key: "_" for key in string.punctuation})
    return x.translate(table)


what_for = 'CL'
convert_cl = makedict()
for file in os.listdir(results_path + what_for + '_P-values/'):
    table_path = results_path + what_for + '_P-values/{}'.format(file)
    table = pd.read_table(table_path)

    if table.empty:
        os.remove(table_path)
        continue

    table.rename({'m_mean_ref': 'es_mean_ref',
                  'm_mean_alt': 'es_mean_alt',
                  'm_mostsig_ref': 'es_mostsig_ref',
                  'm_mostsig_alt': 'es_mostsig_alt', }, axis='columns', inplace=True)
    if what_for == "CL":
        cell_line = convert_cl[file.replace("_common_table.tsv", "")]
        new_table_path = results_path + "FORGTRD/" + what_for + '_P-values/{}.tsv'.format(cell_line)
    else:
        new_table_path = results_path + "FORGTRD/" + what_for + '_P-values/{}'.format(file[:-1 * len("_common_table")])
    print(file, new_table_path)
    table.to_csv(new_table_path, sep="\t", index=False)

for file in os.listdir(results_path + what_for + '_DICTS/'):

    dict_path = results_path + what_for + '_DICTS/{}'.format(file)
    if what_for == "CL":
        new_name = convert_cl[file[:-1*len('_DICT.json')]] + '_DICT.json'
    else:
        new_name = file
    new_dict_path = results_path + "FORGTRD/" + what_for + '_DICTS/{}'.format(new_name)
    print(file, new_name)
    copyfile(dict_path, new_dict_path)
