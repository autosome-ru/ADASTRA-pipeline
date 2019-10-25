#!/bin/bash

ParametersPath="/home/abramov/PARAMETERS/"
ScriptsFolder="/home/abramov/ASB-Project/scripts/"

njobs=$1

parallel --jobs "$njobs" bash "$ScriptsFolder"SNPcalling/ProcessLine.sh :::: "$ParametersPath"/BadVCFs.tsv