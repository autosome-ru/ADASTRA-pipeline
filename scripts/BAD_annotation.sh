#!/bin/bash

source scripts/HELPERS/paths_for_components.py

njobs=$1
remade=$2

if [ "$remade" == --remade ]; then
  adastra make_paths --remade --mode badmaps
  parallel --jobs "$njobs" adastra bad_annotation --remade --base :::: "$parallel_parameters_path"/exp_paths.cfg
else
  adastra make_paths --mode badmaps
  parallel --jobs "$njobs" adastra bad_annotation --base :::: "$parallel_parameters_path"/exp_paths.cfg
fi