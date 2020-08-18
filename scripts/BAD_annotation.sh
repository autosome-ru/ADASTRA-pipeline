#!/bin/bash

source scripts/HELPERS/paths_for_components.py

njobs=$1
python3 -m adastra make_paths --mode badmaps
parallel --jobs "$njobs" python3 -m adastra bad_annotation :::: "$parallel_parameters_path"/exp_paths.cfg