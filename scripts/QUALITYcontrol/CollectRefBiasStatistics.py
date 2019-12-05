import sys
import os.path
import json
import pandas as pd

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import cl_dict_path


name = 'HUES64_embryonic_stem_cells'

with open(cl_dict_path, "r") as read_file:
    cell_lines_dict = json.loads(read_file.readline())
sum_df = None
for align_path in cell_lines_dict[name]:
    if not os.path.isfile(align_path):
        continue
    print(align_path)
    df = pd.read_table(align_path)
    if df.empty:
        print("WHAAAT")
        continue
    if sum_df is None:
        sum_df = df[['#chr', 'pos', 'ref_read_counts', 'alt_read_counts']]
    else:
        sum_df.append(df[['#chr', 'pos', 'ref_read_counts', 'alt_read_counts']])
    print(len(sum_df.index))
print(sum_df['ref_read_counts'].value_counts())
