import pandas as pd
import sys
import os
import statsmodels.stats.multitest

inpDirectory = "/home/abramov/DATA/TF_P-values/"
outDirectory = "/home/abramov/DATA_FOR_MC_FILTER/"


def filterTable(table, mc_tr=10, totc_tr=10, alt=False):
    table = table[table["max_cover"] >= mc_tr]
    table = table[table["total_cover"] >= totc_tr]
    if table.empty:
        return 0
    if alt:
        bool_ar, p_val, _, _ = statsmodels.stats.multitest.multipletests(table["logitp_alt"], alpha=0.05,
                                                                         method='fdr_bh')
        #col = table["logitp_alt"].to_list()
    else:
        bool_ar, p_val, _, _ = statsmodels.stats.multitest.multipletests(table["logitp_ref"], alpha=0.05,
                                                                         method='fdr_bh')
        #col = table["logitp_ref"].to_list()

    return sum(bool_ar)


flag_d = {"totc": "_total_cov", "mc": "_max_cov"}
mc_list = list(range(10, 151, 5))
totc_list = list(range(10, 151, 5))
columns_list = ["#name"] + [str(i) for i in mc_list]
table = pd.DataFrame(columns=columns_list)

alt = bool(int(sys.argv[1]))
if alt:
    alt_str = '_alt'
else:
    alt_str = '_ref'

flag = sys.argv
if flag not in flag_d.keys():
    raise ValueError("{} not in valid modes. {} are valid".format(flag, ",".join(flag_d.keys())))

for filename in os.listdir(inpDirectory):
    with open(inpDirectory + filename, "r") as f:
        noCorrTable = pd.read_table(f)
        if noCorrTable.empty:
            continue
    print("Find statistics for " + filename)
    table.loc[filename, "name"] = filename
    if flag == "mc":
        for MaxCover in mc_list:
            mc_tr = MaxCover
            copy_table = noCorrTable.copy()
            FDR_n = filterTable(copy_table, mc_tr=mc_tr, alt=alt)
            print(FDR_n)
            table.loc[filename, str(mc_tr)] = FDR_n
    if flag == "totc":
        for TotCover in totc_list:
            totc_tr = TotCover
            copy_table = noCorrTable.copy()
            FDR_n = filterTable(copy_table, totc_tr=totc_tr, alt=alt)
            print(FDR_n)
            table.loc[filename, str(totc_tr)] = FDR_n

with open(outDirectory + "statistics_for_TFs" + alt_str + flag_d[flag] + ".tsv", "w") as w:
    table.to_csv(w, sep="\t", index=False)
