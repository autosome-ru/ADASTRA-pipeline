#!/bin/bash

source scripts/HELPERS/paths_for_components.py

njobs=$1

python3 -m adastra make_paths --mode annotation

parallel --jobs "$njobs" bash "$scripts_path"PEAKannotation/ParseMasterLine.sh :::: "$parallel_parameters_path"/exp_paths.cfg