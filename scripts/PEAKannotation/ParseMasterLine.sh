#!/bin/bash

source scripts/HELPERS/soft_configs.cfg
source scripts/HELPERS/paths_for_components.py

PEAKannotationScriptsPath=${scripts_path}"/PEAKannotation/"

IFS=$'\t'
read -ra ADDR <<< "$1"

base_path=${ADDR[0]}
peaks_name=${ADDR[1]}

vcf_path="${base_path}.vcf.gz"
echo "$vcf_path"

echo "Making ${vcf_path}"
echo "Checking exp VCF"
if ! [ -f "$vcf_path" ]; then
  echo "There is no VCF for exp $vcf_path"
  exit 1
fi

if  [ -f "${base_path}.table_annotated.tsv" ]; then
	echo "Remaking annotation for $vcf_path"
else
	echo "Making annotation for $vcf_path first time"
fi

if ! adastra annotate_peaks --base "$base_path"; then
  echo "Failed to annotate $base_path with $peaks_name"
  exit 1
fi

