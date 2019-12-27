import pandas as pd
import sys
import os
import statsmodels.stats.multitest

inpDirectory = "/home/abramov/DATA/TF_P-values/"
outDirectory = "/home/abramov/DATA_FOR_MC_FILTER/"


def filterTable(table, noise_tr=0.0, alt=False):
    table = table[(table["ref_weight_max"] <= noise_tr)]
    if table.empty:
        return 0
    if alt:
        bool_ar, p_val, _, _ = statsmodels.stats.multitest.multipletests(table["logitp_alt"], alpha=0.05,
                                                                         method='fdr_bh')
    else:
        bool_ar, p_val, _, _ = statsmodels.stats.multitest.multipletests(table["logitp_ref"], alpha=0.05,
                                                                         method='fdr_bh')

    return sum(bool_ar)


noise_list = [x/100 for x in range(40)]
columns_list = ["{:.3}".format(i) for i in noise_list]
table = pd.DataFrame(columns=columns_list)
FDRs = {}
for filename in os.listdir(inpDirectory):
    with open(inpDirectory + filename, "r") as f:
        noCorrTable = pd.read_table(f)
        if noCorrTable.empty:
            continue
    print("Find statistics for " + filename)
    for noise in noise_list:
        FDR_n = filterTable(noCorrTable, noise_tr=noise, alt=False)
        try:
            FDRs[noise] += FDR_n
        except KeyError:
            FDRs[noise] = FDR_n
for key in FDRs:
    table.loc["count", "{:.3}".format(key)] = FDRs[key]

with open(outDirectory + "statistics_for_TFs_weight_ref.tsv", "w") as w:
    table.to_csv(w, sep="\t", index=False)
