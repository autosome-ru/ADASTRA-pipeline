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
old_table = pd.read_table(results_path + 'Processed/' + what_for + "_P-values/" + key_name + '_common_table.tsv')

table[['refc_maxdepth', 'altc_maxdepth',
       'BAD_maxdepth', 'm_maxdepth',
       'refc_mostsig', 'altc_mostsig',
       'BAD_mostsig', 'm_mostsig',
       'm_mean_ref', 'm_mean_alt']] = old_table[['refc_maxdepth', 'altc_maxdepth',
                                                 'BAD_maxdepth', 'm_maxdepth',
                                                 'refc_mostsig', 'altc_mostsig',
                                                'BAD_mostsig', 'm_mostsig',
                                                 'm_mean_ref', 'm_mean_alt']]
mc_filter_array = np.array(table['max_cover'] >= 30)
if sum(mc_filter_array) != 0:
    bool_ar_ref, p_val_ref, _, _ = statsmodels.stats.multitest.multipletests(table[mc_filter_array]["logitp_ref"],
                                                                             alpha=0.05, method='fdr_by')
    bool_ar_alt, p_val_alt, _, _ = statsmodels.stats.multitest.multipletests(table[mc_filter_array]["logitp_alt"],
                                                                             alpha=0.05, method='fdr_by')
else:
    p_val_ref = []
    p_val_alt = []
table = table.drop(['fdrp_ref', 'fdrp_alt', 'fdrp_ref_fisher', 'fdrp_alt_fisher'], axis=1)
fdr_by_ref = np.array(['NaN'] * len(table.index))
fdr_by_ref[mc_filter_array] = p_val_ref
table["fdrp_by_ref"] = fdr_by_ref

fdr_by_alt = np.array(['NaN'] * len(table.index))
fdr_by_alt[mc_filter_array] = p_val_alt
table["fdrp_by_alt"] = fdr_by_alt

with open(results_path + 'ProcessedNew/' + what_for + "_P-values/" + key_name + '_common_table.tsv', "w") as w:
    table.to_csv(w, sep="\t", index=False)
