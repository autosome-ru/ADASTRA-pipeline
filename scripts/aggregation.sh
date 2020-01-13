#!/bin/bash

source HELPERS/paths_for_components.py

njobs=$1
indicator=$2

if [ "$indicator" == "--forCL" ];then
  python3 "$scripts_path"PARAMETERS/MakeParametersForAgr.py "CL"
  parallel --jobs "$njobs" python3 "$scripts_path"ASBcalling/Aggregation.py "CL" :::: "$parameters_path"/Agr_parameters.cfg
fi
if [ "$indicator" == "--forTF" ];then
  python3 "$scripts_path"PARAMETERS/MakeParametersForAgr.py "TF"
  parallel --jobs "$njobs" python3 "$scripts_path"ASBcalling/Aggregation.py "TF" :::: "$parameters_path"/Agr_parameters.cfg
fi