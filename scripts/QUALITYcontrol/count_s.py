import pandas as pd
import sys
import os
import statsmodels.stats.multitest

inpDirectory = "/home/abramov/DATA/ProcessedNew/TF_P-values/"
outDirectory = "/home/abramov/DATA_FOR_MC_FILTER/"

counter_001 = 0
counter_005 = 0
for filename in os.listdir(inpDirectory):
    with open(inpDirectory + filename, "r") as f:
        noCorrTable = pd.read_table(f)
        if noCorrTable.empty:
            continue
    print("Find statistics for " + filename)
    counter_005 += len(noCorrTable[(noCorrTable['fdrp_by_ref'] <= 0.05) | (noCorrTable['fdrp_by_alt'] <= 0.05)].index)
    counter_001 += len(noCorrTable[(noCorrTable['fdrp_by_ref'] <= 0.01) | (noCorrTable['fdrp_by_alt'] <= 0.01)].index)
print(counter_005, counter_001)

