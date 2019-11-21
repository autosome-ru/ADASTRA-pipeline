import pandas as pd
import os.path
import sys
import json

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import parameters_path, ploidy_path, ploidy_dict_path
from scripts.HELPERS.helpers import pack

out_path = parameters_path + "cell_lines_BADs.tsv"
actual_ploidy_path = ploidy_path + "Corrected-6/"


def write_BAD(out_buffer, pd_df, max_n, total_n, datasets_n, SNP_n, n_without_SNP, ids):
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
                           datasets_n, SNP_n, n_without_SNP, ",".join(map(str, ids))]))


if __name__ == "__main__":
    with open(ploidy_dict_path, "r") as file:
        cell_lines_dict = json.loads(file.readline())
    with open(out_path, "w") as out:
        out.write(pack(["#cell_line", "mean_BAD_by_SNP", "mean_BAD_by_bp", "bp_max", "bp_total",
                        "number of datasets", "number of SNPs", "number of datasets without SNPs", "geo_encode ID"]))
        sum_table = None
        previous_name = None
        for file_name in sorted(os.listdir(actual_ploidy_path)):
            without_SNP = 0
            split_name = file_name.split("!")
            cell_line_name = split_name[0]
            geo_encode_id = split_name[-1]
            if geo_encode_id == cell_line_name:
                geo_encode_id = "None"
            else:
                geo_encode_id = geo_encode_id.replace("_ploidy.tsv", "")
                print(geo_encode_id)
            if os.stat(ploidy_path + file_name.replace("_ploidy", "")).st_size == 0:
                without_SNP = 1
                cur_SNP_number = 0
            else:
                with open(ploidy_path + file_name.replace("_ploidy", "")) as f:
                    df = pd.read_table(f,  index_col=False)
                    cur_SNP_number = len(df.index)
            cur_l = len([x for x in cell_lines_dict[file_name.split("_ploidy")[0]] if os.path.isfile(x)])
            if without_SNP:
                without_SNP = cur_l
            with open(actual_ploidy_path + file_name) as file:
                table = pd.read_table(file)
                cur_bp_len = (table["end"] - table["start"]).sum()
                if table.empty:
                    without_SNP = cur_l
            if previous_name == cell_line_name:
                bp_len.append(cur_bp_len)
                sum_table = sum_table.append(table)
                geo_encode_list = geo_encode_list.append(geo_encode_id)
                aligns_number += cur_l
                SNP_number += cur_SNP_number
                datasets_without_SNP += without_SNP
            else:
                if sum_table is None:
                    sum_table = table
                    aligns_number = cur_l
                    bp_len = [cur_bp_len]
                    geo_encode_list = [geo_encode_id]
                    previous_name = cell_line_name
                    SNP_number = cur_SNP_number
                    datasets_without_SNP = without_SNP
                else:
                    write_BAD(out, sum_table, max(bp_len), sum(bp_len), aligns_number, SNP_number, datasets_without_SNP,
                              geo_encode_list)
                    datasets_without_SNP = without_SNP
                    sum_table = table
                    bp_len = [cur_bp_len]
                    geo_encode_list = [geo_encode_id]
                    aligns_number = cur_l
                    SNP_number = cur_SNP_number
                    previous_name = cell_line_name
        write_BAD(out, sum_table, max(bp_len), sum(bp_len), aligns_number, SNP_number, datasets_without_SNP,
                  geo_encode_list)
