#!/bin/bash

source HELPERS/paths_for_components.py

njobs=$1

parallel --jobs "$njobs" bash "$scripts_path"PEAKannotation/ParseMasterLine.sh :::: "$parallel_parameters_path"/exp_paths.cfg