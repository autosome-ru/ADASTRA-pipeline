import os
import json
import pandas as pd
from scripts.HELPERS.paths import get_ending
from scripts.HELPERS.paths_for_components import parallel_parameters_path, badmaps_dict_path

out_path = os.path.expanduser('~/badmaps_SNPs.json')


def main():
    with open(badmaps_dict_path, 'r') as read_file:
        d = json.loads(read_file.readline())
    keys = sorted(d.keys())
    out_d = {}
    for key in keys:
        good_count = 0
        total_count = 0
        for value in d[key]:
            path = value + get_ending('annotation')
            if os.path.isfile(path):
                df = pd.read_table(path)
                total_count += len(df.index)
                good_count += len(df[df['ID'] != '.'].index)
        out_d[key] = {'good_count': good_count, 'total_count': total_count}
    with open(out_path, 'w') as file:
        json.dump(out_d, file)
