#!/bin/bash

ParametersListsFolder="/home/abramov/ParallelParameters/stats/"
ScriptsFolder="/home/abramov/ASB-Project/scripts/"

njobs=$1
flag=$2

if [ "$flag" == --merge ]; then
  python3 "$ScriptsFolder"PARAMETERS/MakeParametersForPE.py
  parallel --jobs "$njobs" python3 "$ScriptsFolder"PLOIDYcalling/VCFMerger.py :::: "$ParametersListsFolder"PE_parameters.cfg
fi

if [ "$flag" == --merge ] || [ "$flag" == --ploidy ]; then
  python3 "$ScriptsFolder"PARAMETERS/MakeParametersForPE.py
  python3 "$ScriptsFolder"PARAMETERS/SortParameters.py
	parallel --jobs "$njobs" python3 "$ScriptsFolder"PLOIDYcalling/PloidyEstimation.py :::: "$ParametersListsFolder"PE_parameters.cfg :::: "$ParametersListsFolder"Podgonians.cfg
fi

if [ "$flag" == --merge ] || [ "$flag" == --ploidy ] || [ "$flag" == --aswp ]; then
  python3 "$ScriptsFolder"PARAMETERS/MakeParametersForASWP.py
  parallel --jobs "$njobs" python3 "$ScriptsFolder"CORRELATIONanalysis/Annotate_SNPs_with_ploidy.py :::: "$ParametersListsFolder"ASWP_parameters.cfg
fi

python3 "$ScriptsFolder"PARAMETERS/MakeParametersForCS.py
parallel --jobs "$njobs" python3 "$ScriptsFolder"CORRELATIONanalysis/CorStats.py :::: "$ParametersListsFolder"CS_parameters.cfg

python3 "$ScriptsFolder"CORRELATIONanalysis/JoinThreads.py
