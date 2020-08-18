#!/bin/bash

source scripts/HELPERS/paths_for_components.py
njobs=$1

parallel --delay 80 --jobs "$njobs" \
bash "$scripts_path/SNPcalling/ProcessLine.sh" :::: "${parallel_parameters_path}/sorted_master_list.tsv"
