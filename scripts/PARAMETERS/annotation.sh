#!/bin/bash

ParametersPath="/home/abramov/PARAMETERS/"
ScriptsFolder="/home/abramov/ASB-Project/scripts/"

njobs=$1

parallel --jobs "$njobs" bash "$ScriptsFolder"PEAKannotation/ParseMasterLine.sh :::: "$ParametersPath"/Master-lines.tsv