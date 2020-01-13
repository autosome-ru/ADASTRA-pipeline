#!/bin/bash

source HELPERS/paths_for_components.py

njobs=$1
flag=$2

if [ "$flag" == --merge ]; then
  python3 "$scripts_path"PARAMETERS/MakeParametersForPE.py
  parallel --jobs "$njobs" python3 "$scripts_path"PLOIDYcalling/VCFMerger.py :::: "$parameters_path"PE_parameters.cfg
fi

if [ "$flag" == --merge ] || [ "$flag" == --ploidy ]; then
  python3 "$scripts_path"PARAMETERS/MakeParametersForPE.py
  python3 "$scripts_path"PARAMETERS/SortParameters.py
	parallel --jobs "$njobs" python3 "$scripts_path"PLOIDYcalling/PloidyEstimation.py :::: "$parameters_path"PE_parameters.cfg :::: "$parameters_path"Podgonians.cfg
fi

if [ "$flag" == --merge ] || [ "$flag" == --ploidy ] || [ "$flag" == --aswp ]; then
  python3 "$scripts_path"PARAMETERS/MakeParametersForASWP.py
  parallel --jobs "$njobs" python3 "$scripts_path"CORRELATIONanalysis/Annotate_SNPs_with_ploidy.py :::: "$parameters_path"ASWP_parameters.cfg
fi

python3 "$scripts_path"PARAMETERS/MakeParametersForCS.py
parallel --jobs "$njobs" python3 "$scripts_path"CORRELATIONanalysis/CorStats.py :::: "$parameters_path"CS_parameters.cfg

python3 "$scripts_path"CORRELATIONanalysis/JoinThreads.py
