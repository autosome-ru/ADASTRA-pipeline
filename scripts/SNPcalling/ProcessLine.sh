#!/bin/bash

source ../HELPERS/paths_for_components.py
source ../Configs/CONFIG.cfg

snp_calling_scripts_path=${scripts_path}"/SNPcalling/"
helpers_scripts_path=${scripts_path}"/HELPERS/"

LINE=$1
IFS=$'\t'
read -ra ADDR <<< "$LINE"

exp_name=${ADDR[0]}
align_name=${ADDR[1]}
read_groups=${ADDR[2]}
alignment_download_path=${ADDR[3]}

out_path=${alignments_path}"/$exp_name/"
alignment_full_path=${alignments_path}"/$exp_name/${align_name}.bam"

if [ "$alignment_download_path" == "None" ]; then
    echo "There is no Path for exp ${exp_name}. Checking ${alignments_path}/$exp_name"
    if ! [ -f "$alignment_full_path"]; then
      echo "No data found for ${exp_name}"
      exit 1
    fi
else:
  if ! [ -d ${alignments_path}"$exp_name" ]; then
    if ! mkdir ${alignments_path}"$exp_name"; then
      echo "Failed to make dir $exp_name"
      exit 1
    fi

    echo "Downloading $exp_name"

    if ! bash ${helpers_scripts_path}/DownloadFile.sh "$alignment_download_path" "$alignment_full_path"
    then
      echo "Download failed for $exp_name"
      exit 1
    fi
  fi
fi

echo "Adding read_groups for $exp_name"
if ! bash ${snp_calling_scripts_path}AddReadGroups.sh "$alignment_full_path" "$read_groups"
then
  echo "Failed AddReadGroups $exp_name"
  exit 1
fi

echo "Doing SNPcalling for $exp_name"
if ! bash ${snp_calling_scripts_path}SNPcalling.sh -Exp "$alignment_full_path" -Out "$out_path"
then
  echo "Failed SNPcalling $exp_name"
  exit 1
fi

echo "Cleaning up for ${exp_name}"
rm "$alignment_full_path"
rm "$alignment_full_path.bai"

if [ -f ${out_path}"$align_name.vcf.idx" ]; then
  rm ${out_path}"$align_name.vcf.idx"
fi

if [ -f ${out_path}"$align_name.vcf.gz" ]; then
  rm ${out_path}"$align_name.vcf.gz"
fi

if ! gzip "$out_path$align_name.vcf"; then
	echo "Failed gzip vcf for $exp_name"
	exit 1
fi


