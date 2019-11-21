import pandas as pd
import os.path
import sys
import json

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import parameters_path, ploidy_path, ploidy_dict_path
from scripts.HELPERS.helpers import pack

out_path = parameters_path + "cell_lines_BADs.tsv"
actual_ploidy_path = ploidy_path + "Corrected-6/"


def write_BAD(out_buffer, pd_column, datasets_n, SNP_n, without_SNP):
    out_buffer.write(pack([previous_name, pd_column.mean(), pd_column.median(), datasets_n, SNP_n, without_SNP]))


if __name__ == "__main__":
    with open(ploidy_dict_path, "r") as file:
        cell_lines_dict = json.loads(file.readline())
    with open(out_path, "w") as out:
        out.write(pack(["#cell_line", "mean_BAD", "median_BAD", "number of datasets", "number of SNPs",
                        "number of datasets without SNPs"]))
        sum_table = None
        previous_name = None
        for file_name in sorted(os.listdir(actual_ploidy_path)):
            without_SNP = False
            cell_line_name = file_name.split("!")[0]
            print(ploidy_path + file_name.replace("_ploidy", ""))
            if os.stat(ploidy_path + file_name.replace("_ploidy", "")).st_size == 0:
                without_SNP = True
                cur_SNP_number = 0
            else:
                with open(ploidy_path + file_name.replace("_ploidy", "")) as f:
                    df = pd.read_table(f,  index_col=False)
                    cur_SNP_number = len(df.index)
            cur_l = len([x for x in cell_lines_dict[file_name.split("_ploidy")[0]] if os.path.isfile(x)])
            with open(actual_ploidy_path + file_name) as file:
                table = pd.read_table(file)
            if previous_name == cell_line_name:
                sum_table = sum_table.append(table)
                aligns_number += cur_l
                SNP_number += cur_SNP_number
                datasets_without_SNP += int(without_SNP)
            else:
                if sum_table is None:
                    sum_table = table
                    aligns_number = cur_l
                    previous_name = cell_line_name
                    SNP_number = cur_SNP_number
                    datasets_without_SNP = int(without_SNP)
                else:
                    write_BAD(out, sum_table["BAD"], aligns_number, SNP_number, datasets_without_SNP)
                    datasets_without_SNP = int(without_SNP)
                    sum_table = table
                    aligns_number = cur_l
                    SNP_number = cur_SNP_number
                    previous_name = cell_line_name
        write_BAD(out, sum_table["BAD"], aligns_number, SNP_number, datasets_without_SNP)
