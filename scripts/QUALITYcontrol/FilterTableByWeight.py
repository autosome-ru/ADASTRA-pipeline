import pandas as pd
import sys
import os
import statsmodels.stats.multitest

inpDirectory = "/home/abramov/DATA/TF_P-values/"
outDirectory = "/home/abramov/DATA_FOR_MC_FILTER/"


def filterTable(table, noise_tr=0.0, alt=False):
    table = table[(table["ref_weight_min"] >= noise_tr)]
    if table.empty:
        return 0
    if alt:
        bool_ar, p_val, _, _ = statsmodels.stats.multitest.multipletests(table["logitp_alt"], alpha=0.05,
                                                                         method='fdr_bh')
    else:
        bool_ar, p_val, _, _ = statsmodels.stats.multitest.multipletests(table["logitp_ref"], alpha=0.05,
                                                                         method='fdr_bh')

    return sum(bool_ar)


flag_d = {"totc": "_total_cov", "mc": "_max_cov"}
noise_list = [x/20 for x in range(10)]
columns_list = ["{:.3}".format(i) for i in noise_list]
table = pd.DataFrame(columns=columns_list)

alt = bool(int(sys.argv[1]))
if alt:
    alt_str = '_alt'
else:
    alt_str = '_ref'

for filename in os.listdir(inpDirectory):
    with open(inpDirectory + filename, "r") as f:
        noCorrTable = pd.read_table(f)
        if noCorrTable.empty:
            continue
    print("Find statistics for " + filename)
    for noise in noise_list:
        FDR_n = filterTable(noCorrTable, noise_tr=noise, alt=False)
        table["{:.3}".format(noise)] = FDR_n

with open(outDirectory + "statistics_for_TFs_weight" + alt_str + ".tsv", "w") as w:
    table.to_csv(w, sep="\t")
