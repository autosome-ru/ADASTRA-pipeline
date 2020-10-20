#!/bin/bash

source ../scripts/HELPERS/paths_for_components.py

njobs=$1
adastra make_paths --mode badmaps
parallel --jobs "$njobs" adastra neg_bin_p --base :::: "$parallel_parameters_path"/exp_paths.cfg