from shutil import copy2

import pandas as pd
import json
import os

from scripts.HELPERS.paths_for_components import badmaps_dict_path
from scripts.HELPERS.paths import get_excluded_badmaps_list_path, get_new_badmaps_dict_path, \
    get_badmaps_path_by_validity, get_merged_badmaps_dict_path


def get_bad_dataset_list(log_fdr_tr=5, var_diff_tr=0.05, remake=False):
    filter_df = pd.read_table(get_excluded_badmaps_list_path(remake=remake))
    filter_df = filter_df[
        (filter_df['FDR'] <= 10 ** -log_fdr_tr) &
        (filter_df['dataset_es_var'] - filter_df['ref_es_var'] >= var_diff_tr)
    ]
    return filter_df.apply(lambda row: '{}@{}'.format(row['#Cell_line'], row['Lab']), axis=1).tolist()


def read_json(file):
    with open(file, 'r') as f:
        return json.load(f)


def remake_badmaps_dict(bad_dataset_list):
    """
    :param bad_dataset_list: names of bad datasets <cell_line>@<lab>
    :return: modified badmaps dict and updated old dict
    Select dataset groups from original badmaps dict with poor quality (higher than usual ES variance)
    and split them into separate datasets for subsequent BAD calling at iter 2
    """
    old_dict = read_json(badmaps_dict_path)
    new_dict = {}
    merged_dict = {}
    for dataset, aligns in old_dict.items():
        if dataset in bad_dataset_list and len(aligns) > 1:
            for align in aligns:
                new_dict[dataset.split('@')[0] + '@' + os.path.basename(align)] = [align]
                merged_dict[dataset.split('@')[0] + '@' + os.path.basename(align)] = [align]
        else:
            merged_dict[dataset] = aligns
    return new_dict, merged_dict


def copy_good_badmaps(bad_dataset_list):
    for dataset in os.listdir(os.path.join(get_badmaps_path_by_validity(valid=False), 'CAIC')):
        print(dataset)
        if dataset.split('.')[0] not in bad_dataset_list:
            badmap_path = os.path.join(get_badmaps_path_by_validity(), 'CAIC', dataset)
            new_badmap_path = os.path.join(get_badmaps_path_by_validity(valid=True), 'CAIC', dataset)

            if not os.path.isfile(badmap_path):
                continue

            with open(badmap_path) as src:
                lines = src.readlines()
            if len(lines) <= 1:
                continue
            print(new_badmap_path)
            copy2(badmap_path, new_badmap_path)


def delete_bad_badmaps(bad_dataset_list):
    dir = os.path.join(get_badmaps_path_by_validity(valid=True), 'CAIC')
    for dataset in os.listdir(dir):
        if dataset.split('.')[0] in bad_dataset_list:
            os.remove(os.path.join(dir, dataset))
            continue
        remove = False
        with open(os.path.join(dir, dataset)) as f:
            if len(f.readlines()) <= 1:
                remove = True
        if remove:
            os.remove(os.path.join(dir, dataset))


def main(remake=False):
    bad_dataset_list = get_bad_dataset_list(remake=remake)
    print('Filtered {} datasets'.format(len(bad_dataset_list)))
    if not remake:
        print('iteration 1')
        new_dict, merged_dict = remake_badmaps_dict(bad_dataset_list)
        with open(get_new_badmaps_dict_path(), 'w') as f:
            json.dump(new_dict, f)
        with open(get_merged_badmaps_dict_path(), 'w') as f:
            json.dump(merged_dict, f)
        copy_good_badmaps(bad_dataset_list)
    else:
        print('iteration 2')
        delete_bad_badmaps(bad_dataset_list)
