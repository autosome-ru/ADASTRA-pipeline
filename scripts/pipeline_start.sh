#!/bin/bash

#TODO:
# 1) Dont save parameters files, work with stdin
# 2) Make pipeline start script
# 3) Rename files with "_"
# 4) Create "models_file" for ploidy model and integrate it with corstats

njobs=$1
flag=$2

if [ "$flag" == --BAD ]; then
  bash BAD_annotation.sh "$njobs"
fi

if [ "$flag" == --pvalue ] || [ "$flag" == --BAD ]; then
  bash p_value_count.sh "$njobs"
fi

bash aggregation.sh "$njobs" --forTF

bash aggregation.sh "$njobs" --forCL
