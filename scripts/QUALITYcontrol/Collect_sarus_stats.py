import sys
import os
import numpy as np
import pandas as pd
import json

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths_for_components import parameters_path
from scripts.HELPERS.helpers import states, pack


def get_color(row):
    if abs(row['log_fc']) < np.log10(fc_tr) or abs(row['log_pv']) < -np.log10(fdr_tr):
        return grey_color
    # if row[field + '_ref'] < fdr_tr and row[field + '_alt'] < fdr_tr:
    #     return 'purple'
    if row['log_fc'] * row['log_pv'] > 0:
        return blue_color
    else:
        return red_color


def CollectRedBlue():
    with open(os.path.expanduser("~/PARAMETERS/blue_red_stats.tsv"), "w") as out_file:
        out_file.write(pack(["red", "blue"]))

        for file_name in os.listdir(agr_dir):
            print(file_name)
            pt = pd.read_table(agr_dir + file_name)
            pt = pt.dropna(subset=["motif_fc"])
            if pt.empty:
                continue

            pt = pt[(pt['motif_log_pref'] >= -np.log10(perf_tr)) & (pt['motif_log_palt'] >= -np.log10(perf_tr))]
            pt = pt[~(pt[field + '_alt'].isnull() | pt[field + '_ref'].isnull())]

            pt['log_pv'] = (np.log10(
                pt[[field + '_ref', field + '_alt']]).min(axis=1)) \
                           * np.sign(pt[field + '_alt'] - pt[field + '_ref'])

            pt['log_fc'] = pt['motif_fc']
            pt['col'] = pt.apply(lambda x: get_color(x), axis=1)

            red_counts = len(pt[pt['col'] == red_color].index)
            blue_counts = len(pt[pt['col'] == blue_color].index)
            out_file.write(pack(red_counts, blue_counts))


if __name__ == '__main__':
    agr_dir = os.path.expanduser('~/SARUS/')
    perf_tr = 0.0005
    fc_tr = 4
    fdr_tr = 0.005
    field = 'fdrp_bh'

    blue_color = '#005AB5'  # 1B7837'
    red_color = '#DC3220'  # 762A83'
    grey_color = '#CCCCCC'
    CollectRedBlue()
