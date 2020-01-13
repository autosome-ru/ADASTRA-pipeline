#!/bin/bash

source HELPERS/paths_for_components.py

njobs=$1
flag=$2

parallel --delay 80 --jobs "$njobs" bash "$scripts_path"SNPcalling/ProcessLine.sh "$flag" :::: "$parameters_path"/MasterListForVCFs.tsv
