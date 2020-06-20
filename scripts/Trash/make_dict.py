import json
import sys
import os
import pandas as pd
sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import create_path_from_GTRD_function
from scripts.HELPERS.paths_for_components import parameters_path
from scripts.HELPERS.helpers import remove_punctuation


with open(parameters_path + "CONVERT_CL_NAMES.json") as f:
    dic = json.loads(f.readline())
    old_BAD_dic = {}
    rev_old_BAD_dic = {}
    for key in dic:
        new_val = key.replace("(", "").replace(")", "").replace(" ", "_").replace("/", "_")
        rev_old_BAD_dic[new_val] = dic[key]
        old_BAD_dic[dic[key]] = new_val
    assert len(set(old_BAD_dic.values())) == len(set(old_BAD_dic.keys()))
dirname_new = sys.argv[1]
dirname_old = sys.argv[2]
for filename in os.listdir(dirname_old):
    df_old = pd.read_table(os.path.join(dirname_old, filename))
    filename_new = rev_old_BAD_dic[filename.split("!")[0]] + "!" + filename.split("!")[1]
    df_new = pd.read_table(os.path.join(dirname_new, filename_new))
    if not df_new[['#chr', 'start', 'end', 'BAD']].equals(df_old[['#chr', 'start', 'end', 'BAD']]):
        print('kazino {}'.format(filename_new))
