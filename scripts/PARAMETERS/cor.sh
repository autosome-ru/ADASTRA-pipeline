#!/bin/bash

ParametersListsFolder="/home/abramov/ParallelParameters/"
ScriptsFolder="/home/abramov/ASB-Project/scripts/"

flag=$1

if [ "$flag" == --ploidy ]; then
	bash "$ScriptsFolder"PARAMETERS/aggregation.sh "--forTF" 45  # second parameter - number of jobs
fi

if [ "$flag" == --ploidy ] || [ "$flag" == --aswp ]; then
  python3 "$ScriptsFolder"PARAMETERS/MakeParametersForASWP.py
  python3 "$ScriptsFolder"PARAMETERS/SortParameters.py
  parallel --jobs 40 python3 "$ScriptsFolder"CORRELATIONanalysis/Annotate_SNPs_with_ploidy.py :::: "$ParametersListsFolder"/ASWP_parameters.cfg
fi

python3 "$ScriptsFolder"PARAMETERS/MakeParametersForCS.py
python3 "$ScriptsFolder"PARAMETERS/SortParameters.py
parallel --jobs 40 python3 "$ScriptsFolder"CORRELATIONanalysis/CorStats.py :::: "$ParametersListsFolder"/CS_parameters.cfg

python3 "$ScriptsFolder"CORRELATIONanalysis/JoinThreads.py
