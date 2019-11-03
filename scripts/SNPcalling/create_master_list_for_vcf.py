import sys
import os

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import parameters_path
from scripts.HELPERS.helpers import make_list_for_VCFs

out_path = parameters_path + "MasterListForVCFs_made.tsv"

if __name__ == "__main__":
    make_list_for_VCFs(out_path, condition_function=os.path.isfile)
