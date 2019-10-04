import pandas as pd
import os

inpDirectory = "/home/abramov/RESULTS/TFsForFilter/"
outDirectory = "/home/abramov/RESULTS/TFs_for_PERFECTOS/"
PeakCallersTreshold = -1
MaxCoverTreshold  = 40
SumCoverTreshold = 0


def correctP(p, BonF):
        p = max(p, 0)
        return min(1, p*BonF)

def filterTable(table):
	global PeaksTreshold
	global MaxCoverTreshold
	global SumCoverTreshold
	print("No Bonf table length = {}".format(len(table.index)))
	table = table.loc[table["max_cover"] >= MaxCoverTreshold]
	table = table.loc[table["m_callers"] > PeakCallersTreshold]
	table = table.loc[table["total_cover"] >= SumCoverTreshold]
	BonF=len(table.index)
	table["m_hpref"] = table["m_hpref"].apply(lambda x: correctP(x,BonF))
	table["m_hpalt"] = table["m_hpalt"].apply(lambda x: correctP(x,BonF))
	table["m_fpref"] = table["m_fpref"].apply(lambda x: correctP(x,BonF))
	table["m_fpalt"] = table["m_fpalt"].apply(lambda x: correctP(x,BonF))
	table["m_stpref"] = table["m_stpref"].apply(lambda x: correctP(x,BonF))
	table["m_stpalt"] = table["m_stpalt"].apply(lambda x: correctP(x,BonF))
	table = table.loc[~((table["m_fpref"] == 1) & (table["m_fpalt"] == 1))]
	print("Filtered SNPs {}/{}".format(len(table.index),BonF))
	return table

for filename in os.listdir(inpDirectory):
	with open(inpDirectory + filename,"r") as f:
		noCorrTable = pd.read_table(f)
	f.close()
	print("Filtering table " + filename)
	noCorrTable = filterTable(noCorrTable)
	with open(outDirectory + filename, "w") as w:
		noCorrTable.to_csv(w, sep = "\t", index=False)
	w.close()


