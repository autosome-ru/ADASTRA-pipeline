import pandas as pd
import os.path
import sys
import json

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import parameters_path, ploidy_path, ploidy_dict_path
from scripts.HELPERS.helpers import pack

out_path = parameters_path + "cell_lines_BADs.tsv"
actual_ploidy_path = ploidy_path + "Corrected-6-C-CAIC/"


def write_BAD(out_buffer, pd_df, max_n, total_n, datasets_n, SNP_n, n_without_SNP, IDs):
    print(total_n)
    try:
        mean_by_SNP = (pd_df["BAD"] * pd_df["SNP_count"]).sum()/pd_df["SNP_count"].sum()
    except ZeroDivisionError:
        mean_by_SNP = "nan"
    try:
        mean_by_bp = (pd_df["BAD"] * (pd_df["end"] - pd_df["start"])).sum()/(pd_df["end"] - pd_df["start"]).sum()
    except ZeroDivisionError:
        mean_by_bp = "nan"
    out_buffer.write(pack([previous_name, mean_by_SNP, mean_by_bp, max_n, total_n,
                           datasets_n, SNP_n, n_without_SNP, ",".join(map(str, IDs))]))


if __name__ == "__main__":
    with open(out_path, "w") as out:
        sum_table = None
        previous_name = None
        for file_name in sorted(os.listdir(actual_ploidy_path)):
            out.write(pack(["#segment length in SNPS", "segment length in bp", "file_name"]))
            with open(actual_ploidy_path + file_name) as file:
                for line in file:
                    if line[0] == "#":
                        continue
                    line = line.split()
                    if int(line[7]) > 30000:
                        print("here it comes xD {} ".format(file_name))
                    out.write(pack([line[7], int(line[2])-int(line[1]), file_name]))
