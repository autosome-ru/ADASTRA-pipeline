#!/bin/bash

source HELPERS/paths_for_components.py

njobs=$1

python3 "$scripts_path/"PARAMETERS/make_exp_paths_from_master_list.py 'annotation'

parallel --jobs "$njobs" bash "$scripts_path"PEAKannotation/ParseMasterLine.sh :::: "$parallel_parameters_path"/exp_paths.cfg