import pandas as pd
import sys
import os
import statsmodels.stats.multitest

inpDirectory = "/home/abramov/RESULTS/release-070120/TF_P-values/"
counter_ref = 0
counter_alt = 0
for filename in os.listdir(inpDirectory):
    with open(inpDirectory + filename, "r") as f:
        noCorrTable = pd.read_table(f)
        if noCorrTable.empty:
            continue
    print("Find statistics for " + filename)
    counter_ref += len(noCorrTable[noCorrTable['logitp_ref'] == 0].index)
    counter_alt += len(noCorrTable[noCorrTable['logitp_alt'] == 0].index)

print(counter_alt, counter_ref)
