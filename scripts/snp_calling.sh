#!/bin/bash

source HELPERS/paths_for_components.py
njobs=$1

python3 "$scripts_path/SNPcalling/"sort_columns.py
parallel --delay 80 --jobs "$njobs" \
bash "$scripts_path/SNPcalling/ProcessLine.sh" :::: "${parallel_parameters_path}/sorted_master_list.tsv"
