#!/bin/bash

source scripts/HELPERS/paths_for_components.py

njobs=$1
flag=$2

python3 -m adastra badmaps_params
if [ "$flag" == --merge ]; then
  parallel --jobs "$njobs" python3 -m adastra vcf_merge :::: "$parallel_parameters_path"BE_parameters.cfg
fi

python3 -m adastra sort_params
parallel --jobs "$njobs" python3 -m adastra bad_call :::: "$parallel_parameters_path"BE_parameters.cfg
