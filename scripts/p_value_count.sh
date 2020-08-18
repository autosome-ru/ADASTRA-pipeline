#!/bin/bash

source scripts/HELPERS/paths_for_components.py

njobs=$1
python3 -m adastra make_paths --mode badmaps
parallel --jobs "$njobs" python3 -m adastra neg_bin_p :::: "$parallel_parameters_path"/exp_paths.cfg