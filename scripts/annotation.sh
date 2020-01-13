#!/bin/bash

source HELPERS/paths_for_components.py

njobs=$1

parallel --jobs "$njobs" bash "$scripts_path"PEAKannotation/ParseMasterLine.sh :::: "$parameters_path"/Master-lines.tsv