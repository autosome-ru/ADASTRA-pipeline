#!/bin/bash

source HELPERS/paths_for_components.py

njobs=$1
parallel --jobs "$njobs" python3 "$scripts_path"/ASBcalling/BAD_annotation.py :::: "$parallel_parameters_path"/exp_paths.cfg