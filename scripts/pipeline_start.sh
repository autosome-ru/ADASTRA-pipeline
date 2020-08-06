#!/bin/bash

source HELPERS/paths_for_components.py

njobs=$1
flag=$2

if [ "$flag" == --BAD ]; then
  bash "$scripts_path"/ploidy_est.sh "$njobs" --merge
fi

if [ "$flag" == --BAD ] || [ "$flag" == --NBfit ]; then
  bash "$scripts_path"/BAD_annotation.sh "$njobs"
  python3 "$scripts_path"/FITnoise/CollectRefBiasStatistics.py
  python3 "$scripts_path"/FITnoise/fit_negative_binom_with_weights.py
fi

if [ "$flag" == --pvalue ] || [ "$flag" == --BAD ] || [ "$flag" == --NBfit ]; then
  bash "$scripts_path"/p_value_count.sh "$njobs"
fi
python3 "$scripts_path"/PARAMETERS/MakeDictForAggregation.py TF
python3 "$scripts_path"/PARAMETERS/MakeDictForAggregation.py CL
bash "$scripts_path"/aggregation.sh "$njobs" --forTF
bash "$scripts_path"/aggregation.sh "$njobs" --forCL
