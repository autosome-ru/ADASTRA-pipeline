#!/bin/bash

source ../scripts/HELPERS/paths_for_components.py

njobs=$1
indicator=$2

if [ "$indicator" == "--forCL" ];then
  adastra aggregation_params --for CL
  parallel --jobs "$njobs" adastra aggregation --for CL --name :::: "$parallel_parameters_path"/Agr_parameters.cfg
fi
if [ "$indicator" == "--forTF" ];then
  adastra aggregation_params --for TF
  parallel --jobs "$njobs" adastra aggregation --for TF --name :::: "$parallel_parameters_path"/Agr_parameters.cfg
fi