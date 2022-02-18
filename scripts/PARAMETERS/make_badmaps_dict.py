import requests
import json
import pandas as pd
from scripts.HELPERS.paths import create_path_from_master_list_df
from scripts.HELPERS.paths_for_components import badmaps_dict_path, master_list_path
from scripts.HELPERS.helpers import remove_punctuation, dtype_dict


def find_lab(enc):
    r = requests.get('https://www.encodeproject.org/experiments/{}/?format=json'.format(enc))
    lab = json.loads(r.text)['lab']['@id']
    replicates = json.loads(r.text)['replicates']
    biosamples = sorted(list(set(x['library']['biosample']['@id'] for x in replicates)))
    ret = lab + "_" + ','.join(biosamples)
    return remove_punctuation(ret)


def add_to_dict(d, key, value):
    el = d.get(key, None)
    if el:
        d[key] = el | {value}
    else:
        d[key] = {value}


def add_record(d, row):
    path = create_path_from_master_list_df(row, for_what="base")
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
        # FIXME
        # manually dissemble miss-annotated prostate_cancer
        if key == '22RV1__prostate_carcinoma_@GSE120738':
            key = 'prostate_cancer@{}'.format(row['ALIGNS'])
        add_to_dict(d, key, path)
        return
    elif not pd.isna(row['WG_ENCODE']):
        key = '{}@{}'.format(row['CELLS'], row['WG_ENCODE'])
        add_to_dict(d, key, path)
        return
    raise AssertionError('{} has no ENCODE ID, GEO GSE and wgEncode id'.format(row['#EXP']))


def make_dict(master_list):
    master = pd.read_table(master_list, dtype=dtype_dict)
    master['CELLS'] = master['CELLS'].apply(remove_punctuation)
    d = {}
    df_len = len(master.index)
    for index, row in master.iterrows():
        if df_len > 10 and index % (df_len // 10) == 0:
            print("Made {} Experiments out of {}".format(index, df_len))
        add_record(d, row)
    print("Saving Dictionary")
    for key in d:
        value = d[key]
        sorted_value = sorted(list(value))
        d[key] = sorted_value
    with open(badmaps_dict_path, "w") as write_file:
        json.dump(d, write_file, indent=2)
    print("Dictionary Saved")


def main():
    make_dict(master_list_path)


if __name__ == '__main__':
    main()
