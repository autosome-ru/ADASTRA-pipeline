#!/bin/bash

source scripts/HELPERS/paths_for_components.py

njobs=$1
indicator=$2
remade=$3

if [ "$indicator" == "--forCL" ];then
  if [ "$remade" == --remade ]; then
    adastra aggregation_params --remade --for CL
  else
    adastra aggregation_params --for CL
  fi
  parallel --jobs "$njobs" adastra aggregation --for CL --name :::: "$parallel_parameters_path"/Agr_parameters.cfg
fi
if [ "$indicator" == "--forTF" ];then
  if [ "$remade" == --remade ]; then
    adastra aggregation_params --remade --for TF
  else
    adastra aggregation_params --for TF
  fi
  parallel --jobs "$njobs" adastra aggregation --for TF --name :::: "$parallel_parameters_path"/Agr_parameters.cfg
fi