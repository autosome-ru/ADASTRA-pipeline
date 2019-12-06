import sys
import os.path
import json
import pandas as pd

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import cl_dict_path, parameters_path


# name = [sys.argv[1]]

with open(cl_dict_path, "r") as read_file:
    cell_lines_dict = json.loads(read_file.readline())
out_t = None

for key in cell_lines_dict:
    #if key not in name:
     #   continue
    for align_path in cell_lines_dict[key]:
        if not os.path.isfile(align_path):
            continue
        print(align_path)
        df = pd.read_table(align_path)
        if df.empty:
            continue
        sum_df = df[df['BAD'] == 1][['ref_read_counts', 'alt_read_counts']]

        if out_t is None:
            out_t = pd.DataFrame()
            out_t['ref'] = sum_df['ref_read_counts'].value_counts()
            out_t['alt'] = sum_df['alt_read_counts'].value_counts()
            out_t['allele_reads'] = out_t.index
            out_t = out_t.reset_index(drop=True)
            out_t.fillna(0, inplace=True)
        else:
            tmp_df = pd.DataFrame()
            tmp_df['ref'] = sum_df['ref_read_counts'].value_counts()
            tmp_df['alt'] = sum_df['alt_read_counts'].value_counts()
            tmp_df['allele_reads'] = tmp_df.index
            tmp_df = tmp_df.reset_index(drop=True)
            tmp_df.fillna(0, inplace=True)
            out_t = out_t.append(tmp_df).groupby('allele_reads', as_index=False).sum()

with open(parameters_path + 'bias_statistics.tsv', 'w') as out:
    out_t.to_csv(out, sep="\t")
