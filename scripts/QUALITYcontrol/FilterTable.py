import pandas as pd
import sys
import os
import statsmodels.stats.multitest

inpDirectory = "/home/abramov/DATA/ProcessedNew/TF_P-values/"
outDirectory = "/home/abramov/DATA_FOR_MC_FILTER/"


def filterTable(table, mc_tr=10, totc_tr=10, alt=False):
    table = table[(table["max_cover"] >= mc_tr) | (table["total_cover"] >= totc_tr)]
    if table.empty:
        return 0
    if alt:
        bool_ar, p_val, _, _ = statsmodels.stats.multitest.multipletests(table["logitp_alt"], alpha=0.05,
                                                                         method='fdr_by')
    else:
        bool_ar, p_val, _, _ = statsmodels.stats.multitest.multipletests(table["logitp_ref"], alpha=0.05,
                                                                         method='fdr_by')

    return sum(bool_ar)


def filterTableNAgg(table, mc_tr=10, n=0, alt=False):
    table = table[(table["total_cover"] >= mc_tr) | (table["n_aggregated"] >= n)]
    if table.empty:
        return 0, 0
    if alt:
        bool_ar, p_val, _, _ = statsmodels.stats.multitest.multipletests(table["logitp_alt"], alpha=0.05,
                                                                         method='fdr_bh')
    else:
        bool_ar, p_val, _, _ = statsmodels.stats.multitest.multipletests(table["logitp_ref"], alpha=0.05,
                                                                         method='fdr_bh')

    return sum(bool_ar), len(table.index)

flag_d = {"totc": "_total_cov", "mc": "_max_cov"}
mc_list = list(range(10, 251, 10))
totc_list = list(range(1, 101, 5))
columns_list = [str(i) for i in mc_list]
table = pd.DataFrame(columns=columns_list)

for alt in (True, False):
    if alt:
        alt_str = '_alt'
    else:
        alt_str = '_ref'

    for TotCover in totc_list:
        FDRs = {}
        total_snps_dict = {}
        for filename in os.listdir(inpDirectory):
            with open(inpDirectory + filename, "r") as f:
                noCorrTable = pd.read_table(f)
                if noCorrTable.empty:
                    continue
            print("Find statistics for " + filename)
            for MaxCover in mc_list:
                FDR_n, total_snps = filterTableNAgg(noCorrTable, mc_tr=MaxCover, n=TotCover, alt=alt)
                try:
                    FDRs[MaxCover] += FDR_n
                except KeyError:
                    FDRs[MaxCover] = FDR_n
                try:
                    total_snps_dict[MaxCover] += total_snps
                except KeyError:
                    total_snps_dict[MaxCover] = total_snps
        for MaxCover in mc_list:
            table.loc[str(TotCover), str(MaxCover)] = str(FDRs[MaxCover]) + ',' + str(total_snps_dict[MaxCover])

    with open(outDirectory + "statistics_for_TFs" + alt_str + ".tsv", "w") as w:
        table.to_csv(w, sep="\t")
