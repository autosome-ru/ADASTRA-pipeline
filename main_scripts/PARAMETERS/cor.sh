#!/bin/bash

ParametersListsFolder="/home/abramov/ParallelParameters/"
ScriptsFolder="/home/abramov/ASB-Project/main_scripts/"

flag=$1

if [ "$flag" == --ploidy ]; then
	parallel --jobs 80 python3 "$ScriptsFolder"PLOIDYcalling/PloidyEstimation.py :::: "$ParametersListsFolder"/PE_parameters.cfg
fi

python3 "$ScriptsFolder"PARAMETERS/MakeParametersForASWP.py

parallel --jobs 40 python3 "$ScriptsFolder"CORRELATIONanalysis/Annotate_SNPs_with_ploidy.py :::: "$ParametersListsFolder"/ASWP_parameters.cfg

python3 "$ScriptsFolder"CORRELATIONanalysis/CorStats.py