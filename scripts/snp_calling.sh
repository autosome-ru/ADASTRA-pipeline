#!/bin/bash

source HELPERS/paths_for_components.py
njobs=$1
flag=$2

snp_calling_path="$scripts_path"/SNPcalling/
bash "$snp_calling_path"CreateReference.sh --RefFolder "$reference_path" --RefGenome "$genome_path"
python3 "$snp_calling_path"create_master_list_for_vcf.py

parallel --delay 80 --jobs "$njobs" bash "$snp_calling_path"ProcessLine.sh "$flag" :::: "$configs_path"/MasterListForVCFs.tsv
