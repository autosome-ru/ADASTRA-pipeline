import pandas as pd
import os.path
import sys

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import parameters_path, ploidy_path
from scripts.HELPERS.helpers import pack

out_path = parameters_path + "cell_lines_BADs.tsv"
actual_ploidy_path = ploidy_path + "Corrected-6/"

if __name__ == "__main__":
    with open(out_path, "w") as out:
        out.write(pack(["#cell_line", "mean_BAD", "median_BAD"]))
        sum_table = None
        previous_name = None
        for file_name in sorted(os.listdir(actual_ploidy_path)):
            cell_line_name = file_name.split("!")[0]
            print(cell_line_name)
            with open(actual_ploidy_path + file_name) as file:
                table = pd.read_table(file)
            if previous_name == cell_line_name:
                sum_table = sum_table.append(table)
            else:
                if sum_table is None:
                    sum_table = table
                    previous_name = cell_line_name
                else:
                    mean_BAD = sum_table["BAD"].mean()
                    median_BAD = sum_table["BAD"].median()
                    out.write(pack([previous_name, mean_BAD, median_BAD]))
                    sum_table = table
                    previous_name = cell_line_name
        mean_BAD = sum_table["BAD"].mean()
        median_BAD = sum_table["BAD"].median()
        out.write(pack([previous_name, mean_BAD, median_BAD]))
