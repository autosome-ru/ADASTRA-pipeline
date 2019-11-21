import pandas as pd
import os.path
import sys
import json

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import parameters_path, ploidy_path, cl_dict_path
from scripts.HELPERS.helpers import pack

out_path = parameters_path + "cell_lines_BADs.tsv"
actual_ploidy_path = ploidy_path + "Corrected-6/"


def write_BAD(out_buffer, pd_column, datasets_n):
    out_buffer.write(pack([previous_name, pd_column.mean(), pd_column.median(), datasets_n, len(pd_column.index)]))


if __name__ == "__main__":
    with open(cl_dict_path, "r") as file:
        cell_lines_dict = json.loads(file.readline())
    with open(out_path, "w") as out:
        out.write(pack(["#cell_line", "mean_BAD", "median_BAD", "number of datasers", "number of SNPs"]))
        sum_table = None
        previous_name = None
        for file_name in sorted(os.listdir(actual_ploidy_path)):
            cell_line_name = file_name.split("!")[0]
            try:
                cur_l = len([x for x in cell_lines_dict[cell_line_name] if os.path.isfile(x)])
            except KeyError:
                cur_l = "Bug"
            with open(actual_ploidy_path + file_name) as file:
                table = pd.read_table(file)
            if previous_name == cell_line_name:
                sum_table = sum_table.append(table)
                aligns_number += cur_l
            else:
                if sum_table is None:
                    sum_table = table
                    aligns_number = cur_l
                    previous_name = cell_line_name
                else:
                    write_BAD(out, sum_table["BAD"], aligns_number)
                    sum_table = table
                    aligns_number = cur_l
                    previous_name = cell_line_name
        write_BAD(out, sum_table["BAD"], aligns_number)
