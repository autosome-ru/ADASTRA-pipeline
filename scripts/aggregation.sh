#!/bin/bash

source HELPERS/paths_for_components.py

njobs=$1
indicator=$2

if [ "$indicator" == "--forCL" ];then
  python3 -m ADASTRA aggregation_params CL
  parallel --jobs "$njobs" python3 -m ADASTRA aggregation CL :::: "$parallel_parameters_path"/Agr_parameters.cfg
fi
if [ "$indicator" == "--forTF" ];then
    python3 -m ADASTRA aggregation_params TF
  parallel --jobs "$njobs" python3 -m ADASTRA aggregation TF :::: "$parallel_parameters_path"/Agr_parameters.cfg
fi