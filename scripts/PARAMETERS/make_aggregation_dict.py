import json
import os

import pandas as pd

from scripts.HELPERS.helpers import check_if_in_expected_args, remove_punctuation, dtype_dict
from scripts.HELPERS.paths import create_path_from_master_list_df
from scripts.HELPERS.paths_for_components import configs_path, master_list_path


def makedict(what_for):
    d = dict()
    check_if_in_expected_args(what_for)
    master_df = pd.read_table(master_list_path, dtype=dtype_dict)
    master_df['path'] = master_df.apply(lambda x: create_path_from_master_list_df(x, for_what='p-value'))

    for index, row in master_df.iterrows():
        if what_for == "TF":
            try:
                d[row['TF_UNIPROT_ID']].append(row['path'])
            except KeyError:
                d[row['TF_UNIPROT_ID']] = [row['path']]
        if what_for == "CL":
            cell_line = remove_punctuation(row['CELLS'])
            try:
                d[cell_line].append(row['path'])
            except KeyError:
                d[cell_line] = [row['path']]
    with open(os.path.join(configs_path, what_for + "_DICT.json"), "w") as write_file:
        json.dump(d, write_file)
    print("Dictionary Saved")


def main():
    for ind in 'TF', 'CL':
        makedict(ind)


if __name__ == '__main__':
    main()
