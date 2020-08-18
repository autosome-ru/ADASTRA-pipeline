#!/bin/bash

source scripts/HELPERS/paths_for_components.py

njobs=$1
indicator=$2

if [ "$indicator" == "--forCL" ];then
  adastra aggregation_params CL
  parallel --jobs "$njobs" adastra aggregation CL :::: "$parallel_parameters_path"/Agr_parameters.cfg
fi
if [ "$indicator" == "--forTF" ];then
  adastra aggregation_params TF
  parallel --jobs "$njobs" adastra aggregation TF :::: "$parallel_parameters_path"/Agr_parameters.cfg
fi