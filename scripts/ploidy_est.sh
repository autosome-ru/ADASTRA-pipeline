#!/bin/bash

source HELPERS/paths_for_components.py

njobs=$1
flag=$2

python3 "$scripts_path"PARAMETERS/MakeParametersForPE.py
if [ "$flag" == --merge ]; then
  parallel --jobs "$njobs" python3 "$scripts_path"PLOIDYcalling/VCFMerger.py :::: "$parallel_parameters_path"PE_parameters.cfg
fi

python3 "$scripts_path"PARAMETERS/SortParameters.py
parallel --jobs "$njobs" python3 "$scripts_path"PLOIDYcalling/PloidyEstimation.py :::: "$parallel_parameters_path"PE_parameters.cfg :::: "$parallel_parameters_path"Podgonians.cfg
