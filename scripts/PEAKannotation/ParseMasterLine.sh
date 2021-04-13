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
for peak in $intervals_path/*; do
  peak_type="$(basename "$(dirname "$peak")")"
  peak_path="$peak/$peaks_name.interval"
  if [ -f "$peak_path" ]; then
    if ! adastra check_pos_peaks --peak "$peak_path" --out "${base_path}.${peak_type}.bed" --type $peak_type; then
      echo  "Failed check pos peaks"
      exit 1
    fi

    if ! bedtools sort -i "${base_path}.$peak_type.bed" > "${base_path}.${peak_type}.bed.sorted"
      then
        echo "Failed to sort $peak_type peaks"
        exit 1
    fi
    rm "${base_path}.$peak_type.bed"
  fi
done
if ! adastra annotate_peaks --base "$base_path"; then
  echo "Failed to annotate $base_path with $peaks_name"
  exit 1
fi

for peak in $intervals_path; do
  peak_type="$(basename "$(dirname "$peak")")"
  if [ -f "${base_path}.$peak_type.bed.sorted" ]; then
    rm "${base_path}.$peak_type.bed.sorted"
  fi
done
