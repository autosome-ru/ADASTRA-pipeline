#!/bin/bash

script_path="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
source $script_path/HELPERS/paths_for_components.py

njobs=$1

ls -1 ${results_path}/TF_P-values | sed -e 's/\.tsv$//' | parallel --jobs "$njobs" bash "$scripts_path"/SARUSannotation/start_sarus.sh