import sys
import os.path
import json
import pandas as pd

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import cl_dict_path, parameters_path


name = sys.argv[1]

with open(cl_dict_path, "r") as read_file:
    cell_lines_dict = json.loads(read_file.readline())
sum_df = None
for align_path in cell_lines_dict[name]:
    if not os.path.isfile(align_path):
        continue
    print(align_path)
    df = pd.read_table(align_path)
    if df.empty:
        continue
    if sum_df is None:
        sum_df = df[['ref_read_counts', 'alt_read_counts']]
    else:
        sum_df = sum_df.append(df[['ref_read_counts', 'alt_read_counts']])
with open(parameters_path + 'bias_statistics.tsv', 'w') as out:
    out_t = pd.DataFrame()
    ser = sum_df['ref_read_counts'].value_counts()
    out_t['ref'] = sum_df['ref_read_counts'].value_counts()
    out_t['alt'] = sum_df['alt_read_counts'].value_counts()
    out_t['allele_reads'] = out_t.index
    out_t.reset_index(drop=True)
    out_t.fillna(0, inplace=True)
    out_t.to_csv(out, sep="\t")
