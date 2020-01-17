#!/bin/bash

source HELPERS/paths_for_components.py

njobs=$1

parallel --jobs "$njobs" bash "$scripts_path"SarusFC.sh ::: "${results_path}TF_P-values/"*
