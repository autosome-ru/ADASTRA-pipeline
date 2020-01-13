#!/bin/bash

source HELPERS/paths_for_components.py

njobs=$1
python3 "$scripts_path"PARAMETERS/MakeParametersForPvC.py
parallel --jobs "$njobs" python3 "$scripts_path"ASBcalling/BAD_annotation.py :::: "$parameters_path"/PvC_parameters.cfg