#!/bin/bash

source scripts/HELPERS/paths_for_components.py

njobs=$1
flag=$2

adastra badmaps_params
if [ "$flag" == --merge ]; then
  parallel --jobs "$njobs" adastra vcf_merge :::: "$parallel_parameters_path"/BE_parameters.cfg
fi

adastra sort_params
parallel --jobs "$njobs" adastra bad_call :::: "$parallel_parameters_path"/BE_parameters.cfg
