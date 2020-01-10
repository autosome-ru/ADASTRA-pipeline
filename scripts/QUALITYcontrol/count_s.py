import pandas as pd
import sys
import os
import statsmodels.stats.multitest

inpDirectory = "/home/abramov/RESULTS/release-100120_Dipper/TF_P-values/"
outDirectory = "/home/abramov/DATA_FOR_MC_FILTER/"

counter_001 = 0
counter_005 = 0
datasets_counter = 0
datasets_counter_fdr = 0
all_snps_counter = 0
for filename in os.listdir(inpDirectory):
    with open(inpDirectory + filename, "r") as f:
        noCorrTable = pd.read_table(f)
        if noCorrTable.empty:
            continue
    datasets_counter += 1
    print("Find statistics for " + filename)

    all_snps_counter += len(noCorrTable.index)
    cur_counter_005 = len(noCorrTable[(noCorrTable['fdrp_by_ref'] <= 0.05)].index) + len(noCorrTable[noCorrTable['fdrp_by_alt'] <= 0.05].index)
    cur_counter_001 = len(noCorrTable[(noCorrTable['fdrp_by_ref'] <= 0.01)].index) + len(noCorrTable[noCorrTable['fdrp_by_alt'] <= 0.01].index)
    counter_001 += cur_counter_001
    counter_005 += cur_counter_005
    if cur_counter_005 > 0:
        datasets_counter_fdr += 1


print('0.05 FDR BY snp: {}\n'
      '0.01 FDR BY snp: {}\n'
      'all_snps: {}\n'
      'datasets: {}\n'
      'datasets with asb snps: {}\n'.format(counter_005, counter_001, all_snps_counter, datasets_counter, datasets_counter_fdr))

