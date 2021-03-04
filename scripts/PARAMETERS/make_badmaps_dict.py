import requests
import json
import pandas as pd
import os
from scripts.HELPERS.paths import create_path_from_master_list_df
from scripts.HELPERS.paths_for_components import badmaps_dict_path, master_list_path
from scripts.HELPERS.helpers import remove_punctuation, dtype_dict


def find_lab(enc):
    r = requests.get('https://www.encodeproject.org/experiments/{}/?format=json'.format(enc))
    lab = json.loads(r.text)['lab']['@id']
    biosample = json.loads(r.text)['replicates'][0]['library']['biosample']['@id']
    ret = lab + "_" + biosample
    return remove_punctuation(ret)


def add_to_dict(d, key, value):
    el = d.get(key, None)
    if el:
        d[key] = el | {value}
    else:
        d[key] = {value}


def add_record(d, row, is_chip=False):
    path = create_path_from_master_list_df(row, for_what="base", is_chip=is_chip)
    if not pd.isna(row['ENCODE']):
        Lab = find_lab(row['ENCODE'])
        if Lab:
            key = '{}@{}'.format(row['CELLS'], Lab)
            add_to_dict(d, key, path)
            return
        else:
            raise AssertionError('Lab not found')
    elif not pd.isna(row['GEO']):
        key = '{}@{}'.format(row['CELLS'], row['GEO'])
        add_to_dict(d, key, path)
        return
    elif not pd.isna(row['WG_ENCODE']):
        key = '{}@{}'.format(row['CELLS'], row['WG_ENCODE'])
        add_to_dict(d, key, path)
        return
    raise AssertionError('{} has no ENCODE ID, GEO GSE and wgEncode id'.format(row['#EXP']))


def make_dict(master_list_dnase, master_list_chip):
    master_chip = None
    if master_list_chip:
        master_chip = pd.read_table(master_list_chip, dtype=dtype_dict)
        master_chip['CELLS'] = master_chip['CELLS'].apply(remove_punctuation)
    master_dnase = pd.read_table(master_list_dnase, dtype=dtype_dict)
    master_dnase['CELLS'] = master_dnase['CELLS'].apply(remove_punctuation)
    d = {}
    if master_chip:
        chip_len = len(master_chip.index)
        for index, row in master_chip.iterrows():
            if chip_len > 10 and index % (chip_len // 10) == 0:
                print("Made {} Experiments out of {}".format(index, chip_len))
            add_record(d, row, is_chip=True)
    dnase_len = len(master_dnase.index)
    for index, row in master_dnase.iterrows():
        if dnase_len > 10 and index % (dnase_len // 10) == 0:
            print("Made {} Experiments out of {}".format(index, dnase_len))
        add_record(d, row, is_chip=False)

    print("Saving Dictionary")
    for key in d:
        value = d[key]
        sorted_value = sorted(list(value))
        d[key] = sorted_value
    with open(badmaps_dict_path, "w") as write_file:
        json.dump(d, write_file)
    print("Dictionary Saved")


def main():
    make_dict(master_list_path, os.path.expanduser('~/Configs/master-chip.txt'))


if __name__ == '__main__':
    main()
