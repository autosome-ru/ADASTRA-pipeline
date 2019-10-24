#!/bin/bash

ParametersListsFolder="/home/abramov/ParallelParameters/"
ScriptsFolder="/home/abramov/ASB-Project/scripts/"

njobs=$1
flag=$2

python3 "$ScriptsFolder"PARAMETERS/MakeParametersForPE.py
if [ "$flag" == --merge ]; then
  parallel --jobs "$njobs" python3 "$ScriptsFolder"PLOIDYcalling/VCFMerger.py :::: "$ParametersListsFolder"PE_parameters.cfg
fi

python3 "$ScriptsFolder"PARAMETERS/SortParameters.py
parallel --jobs "$njobs" python3 "$ScriptsFolder"PLOIDYcalling/PloidyEstimation.py :::: "$ParametersListsFolder"/PE_parameters.cfg