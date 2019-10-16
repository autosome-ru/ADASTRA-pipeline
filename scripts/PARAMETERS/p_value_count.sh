#!/bin/bash

ParametersListsFolder="/home/abramov/ParallelParameters/"
ScriptsFolder="/home/abramov/ASB-Project/main_scripts/"
njobs=$1
python3 "$ScriptsFolder"PARAMETERS/MakeParametersForPvC.py
parallel --jobs "$njobs" python3 "$ScriptsFolder"ASBcalling/Pcounter.py :::: "$ParametersListsFolder"/PvC_parameters.cfg