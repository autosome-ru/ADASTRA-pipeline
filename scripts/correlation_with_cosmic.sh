#!/bin/bash

source scripts/HELPERS/paths_for_components.py

njobs=$1
flag=$2

if [ "$flag" == --bad-call ]; then
  adastra badmaps_params
  parallel --jobs "$njobs" adastra vcf_merge --group :::: "$parallel_parameters_path"/BE_parameters.cfg
  adastra sort_params
  parallel --jobs "$njobs" adastra bad_call --group :::: "$parallel_parameters_path"/BE_parameters.cfg
fi


if [ "$flag" == --bad-call ] || [ "$flag" == --annotate ]; then
  adastra annotation_params
  parallel --jobs "$njobs" adastra annotate_snps_for_correlation --base :::: "$parallel_parameters_path"/ASWP_parameters.cfg
fi

if [ "$flag" == --bad-call ] || [ "$flag" == --annotate ] || [ "$flag" == --correlation ]; then

  adastra correlation_params
  parallel --jobs "$njobs" adastra cosmic_correlation --base :::: "$parallel_parameters_path"/CS_parameters.cfg

  adastra join_correlation_threads
fi
