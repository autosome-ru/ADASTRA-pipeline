#!/bin/bash

source HELPERS/paths_for_components.py

#TODO:
# 1) Dont save parameters files, work with stdin
# 2) Make pipeline start script
# 3) Rename files with "_"
# 4) Create "models_file" for ploidy model and integrate it with corstats

njobs=$1
flag=$2

if [ "$flag" == --BAD ]; then
  bash ploidy_est.sh "$njobs" --merge
  bash BAD_annotation.sh "$njobs"
fi

if [ "$flag" == --BAD ] || [ "$flag" == --NBfit ]; then
  python3 "$scripts_path"QUALITYcontrol/CollectRefBiasStatistics.py
  python3 "$scripts_path"FITnoise/fit_negative_binom_with_weights.py
fi

if [ "$flag" == --pvalue ] || [ "$flag" == --BAD ] || [ "$flag" == --NBfit ]; then
  bash p_value_count.sh "$njobs"
fi
python3 "$scripts_path"PARAMETERS/MakeDictForAggregation.py TF
python3 "$scripts_path"PARAMETERS/MakeDictForAggregation.py CL
bash aggregation.sh "$njobs" --forTF

bash aggregation.sh "$njobs" --forCL
