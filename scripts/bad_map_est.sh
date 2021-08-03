#!/bin/bash

source scripts/HELPERS/paths_for_components.py

njobs=$1
flag=$2
redo=$3

if [ "$redo" == --remake ]; then
  adastra badmaps_params --remake
else
  adastra badmaps_params
fi
if [ "$flag" == --merge ]; then
  if [ "$redo" == --remake ]; then
    parallel --jobs "$njobs" adastra vcf_merge --remake --group :::: "$parallel_parameters_path"/BE_parameters.cfg
  else
    parallel --jobs "$njobs" adastra vcf_merge --group :::: "$parallel_parameters_path"/BE_parameters.cfg
  fi
fi

adastra sort_params
if [ "$redo" == --remake ]; then
  parallel --jobs "$njobs" adastra bad_call --remake --group :::: "$parallel_parameters_path"/BE_parameters.cfg
else
  parallel --jobs "$njobs" adastra bad_call --group :::: "$parallel_parameters_path"/BE_parameters.cfg
fi