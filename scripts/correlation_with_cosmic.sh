#!/bin/bash

source scripts/HELPERS/paths_for_components.py

njobs=$1
flag=$2
redo=$3

if [ "$flag" == --bad-call ]; then
  adastra badmaps_params
  if [ "$redo" == --remake ]; then
    parallel --jobs "$njobs" adastra vcf_merge --remake --group :::: "$parallel_parameters_path"/BE_parameters.cfg
  else
    parallel --jobs "$njobs" adastra vcf_merge --group :::: "$parallel_parameters_path"/BE_parameters.cfg
  fi
  adastra sort_params
  if [ "$redo" == --remake ]; then
    parallel --jobs "$njobs" adastra bad_call --remake --group :::: "$parallel_parameters_path"/BE_parameters.cfg
  else
    parallel --jobs "$njobs" adastra bad_call --group :::: "$parallel_parameters_path"/BE_parameters.cfg
  fi
fi


if [ "$flag" == --bad-call ] || [ "$flag" == --annotate ]; then
  adastra annotation_params
  if [ "$redo" == --remake ]; then
    parallel --jobs "$njobs" adastra annotate_snps_for_correlation --remake --base :::: "$parallel_parameters_path"/ASWP_parameters.cfg
  else
    parallel --jobs "$njobs" adastra annotate_snps_for_correlation --base :::: "$parallel_parameters_path"/ASWP_parameters.cfg
  fi
fi

if [ "$flag" == --bad-call ] || [ "$flag" == --annotate ] || [ "$flag" == --correlation ]; then
  if [ "$redo" == --remake ]; then
    adastra correlation_params --remake
    parallel --jobs "$njobs" adastra cosmic_correlation --remake --base :::: "$parallel_parameters_path"/CS_parameters.cfg
    adastra join_correlation_threads --remake
  else
    adastra correlation_params
    parallel --jobs "$njobs" adastra cosmic_correlation --base :::: "$parallel_parameters_path"/CS_parameters.cfg
    adastra join_correlation_threads
  fi
fi
