import sys
import os.path
from statistics import median_grouped
from scipy import stats
import numpy as np
import json
import statsmodels.stats.multitest
import pandas as pd

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import results_path, cl_dict_path, tf_dict_path, parameters_path
from scripts.HELPERS.helpers import callers_names, unpack, pack, check_if_in_expected_args, expected_args, states

what_for = sys.argv[1]  # "TF" or "CL" arguments are expected
check_if_in_expected_args(what_for)
key_name = sys.argv[2]
print("Counting FDR for {} {}".format(what_for, key_name))
table = pd.read_table(results_path + what_for + "_P-values/" + key_name + '_common_table.tsv')
if table.empty:
    sys.exit(0)

bool_ar_ref, p_val_ref, _, _ = statsmodels.stats.multitest.multipletests(table["logitp_ref"],
                                                                         alpha=0.05, method='fdr_by')
bool_ar_alt, p_val_alt, _, _ = statsmodels.stats.multitest.multipletests(table["logitp_alt"],
                                                                         alpha=0.05, method='fdr_by')
table["fdrp_by_ref"] = pd.Series(p_val_ref)
table["fdrp_by_alt"] = pd.Series(p_val_alt)


with open(results_path + 'Processed/' + what_for + "_P-values/" + key_name + '_common_table.tsv', "w") as w:
    table.to_csv(w, sep="\t", index=False)
