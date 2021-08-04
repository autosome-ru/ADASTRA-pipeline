#!/bin/bash

source scripts/HELPERS/paths_for_components.py

njobs=$1
remade=$2
adastra make_paths --mode badmaps
if [ "$remade" == --remade ]; then
  parallel --jobs "$njobs" adastra bad_annotation --remade --base :::: "$parallel_parameters_path"/exp_paths.cfg
else
  parallel --jobs "$njobs" adastra bad_annotation --base :::: "$parallel_parameters_path"/exp_paths.cfg
fi