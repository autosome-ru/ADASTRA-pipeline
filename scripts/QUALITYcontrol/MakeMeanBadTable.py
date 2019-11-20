import pandas as pd
import os.path
import sys

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import parameters_path, ploidy_path
from scripts.HELPERS.helpers import pack

out_path = parameters_path + "cell_lines_BADs.tsv"

if __name__ == "__main__":
    with open(out_path, "w") as out:
        sum_table = None
        previous_name = None
        for file_name in sorted(os.listdir(ploidy_path)):
            cell_line_name = file_name.split("!")[0]
            if previous_name == cell_line_name:
                table = pd.read_table(ploidy_path + file_name)
                sum_table = sum_table.append(table)
            else:
                if sum_table is None:
                    sum_table = pd.read_table(ploidy_path + file_name)
                else:
                    mean_BAD = sum_table["BAD"].mean()
                    median_BAD = sum_table["BAD"].median()
                    out.write(pack([previous_name, mean_BAD, median_BAD]))
                    sum_table = pd.read_table(ploidy_path + file_name)
                    previous_name = cell_line_name
