#!/bin/bash

source scripts/HELPERS/paths_for_components.py

njobs=$1

ls -1 ${results_path}/TF_P-values | sed -e 's/\.tsv$//' | parallel --jobs "$njobs" bash "$scripts_path"/SARUSannotation/start_sarus.sh