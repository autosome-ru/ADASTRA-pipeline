#!/bin/bash

ParametersListsFolder="/home/abramov/ParallelParameters/"
ScriptsFolder="/home/abramov/ASB-Project/main_scripts/"
indicator=$1
njobs=$2

if [ $indicator == "--forCL" ];then
  python3 "$ScriptsFolder"PARAMETERS/MakeParametersForAgr.py "CL"
  parallel --jobs "$njobs" python3 "$ScriptsFolder"ASBcalling/NoBonfSoftAgregation.py "CL" :::: "$ParametersListsFolder"/Agr_parameters.cfg
fi
if [ $indicator == "--forTF" ];then
  python3 "$ScriptsFolder"PARAMETERS/MakeParametersForAgr.py "TF"
  parallel --jobs "$njobs" python3 "$ScriptsFolder"ASBcalling/NoBonfSoftAgregation.py "TF" :::: "$ParametersListsFolder"/Agr_parameters.cfg
fi