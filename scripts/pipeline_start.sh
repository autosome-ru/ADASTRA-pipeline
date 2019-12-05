#!/bin/bash

#TODO:
# 1) Dont save parameters files, work with stdin
# 2) Make pipeline start script
# 3) Rename files with "_"
# 4) Create "models_file" for ploidy model and integrate it with corstats
# 5) move bash files from PARAMETERS dir

ScriptsFolder="/home/abramov/ASB-Project/scripts/"

njobs=$1
flag=$2

python3 "$ScriptsFolder"PARAMETERS/MakeParametersForPE.py
if [ "$flag" == --pvalue ]; then
  bash "$ScriptsFolder"PARAMETERS/p_value_count.sh "$njobs"
fi

bash "$ScriptsFolder"PARAMETERS/aggregation.sh "$njobs" --forTF

bash "$ScriptsFolder"PARAMETERS/aggregation.sh "$njobs" --forCL
