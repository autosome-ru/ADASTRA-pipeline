#!/bin/bash

DownloadPath=$1
BamPath=$2

if ! scp -T -P 1300 autosome@localhost:"$download_path" "$dir_path"
then
  echo "Failed to download $download_path in $dir_path"
  exit 1
fi

python3 ${scripts_path}/HELPERS/make_bam_list "$dir_path" $download_path

if ! samtools merge -b "$dir_path/bam_list.txt" "$dir_path"/${bam_name}
then
  echo "Failed to merge vcfs for $download_path"
  exit 1
fi

