#!/bin/bash

ParametersListsFolder="/home/abramov/ParallelParameters/"
ScriptsFolder="/home/abramov/ASB-Project/main_scripts/"
njobs=$1
python3 "$ScriptsFolder"PARAMETERS/MakeParametersForPE.py
parallel --jobs "$njobs" python3 "$ScriptsFolder"PLOIDYcalling/PloidyEstimation.py :::: "$ParametersListsFolder"/PE_parameters.cfg