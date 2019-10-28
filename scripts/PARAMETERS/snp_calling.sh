#!/bin/bash

ParametersPath="/home/abramov/PARAMETERS/"
ScriptsFolder="/home/abramov/ASB-Project/scripts/"

njobs=$1
flag=$2

parallel --tmux --delay 80 --jobs "$njobs" bash "$ScriptsFolder"SNPcalling/ProcessLine.sh "$flag" :::: "$ParametersPath"/BadVCF.tsv

