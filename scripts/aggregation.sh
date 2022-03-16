#!/bin/bash

source scripts/HELPERS/paths_for_components.py

njobs=$1
indicator=$2
remade=$3

adastra aggregation_dict

if [ "$indicator" == "--forCL" ];then
  if [ "$remade" == "--remade" ]; then
    adastra aggregation_params --remade --for CL
    parallel --jobs "$njobs" adastra aggregation --remade --for CL --name :::: "$parallel_parameters_path"/Agr_parameters.cfg
  else
    adastra aggregation_params --for CL
    parallel --jobs "$njobs" adastra aggregation --for CL --name :::: "$parallel_parameters_path"/Agr_parameters.cfg
  fi
fi
if [ "$indicator" == "--forTF" ];then
  if [ "$remade" == "--remade" ]; then
    adastra aggregation_params --remade --for TF
    parallel --jobs "$njobs" adastra aggregation --remade --for TF --name :::: "$parallel_parameters_path"/Agr_parameters.cfg
  else
    adastra aggregation_params --for TF
    parallel --jobs "$njobs" adastra aggregation --for TF --name :::: "$parallel_parameters_path"/Agr_parameters.cfg
  fi
fi