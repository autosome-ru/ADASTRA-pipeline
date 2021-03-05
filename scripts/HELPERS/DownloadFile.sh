#!/bin/bash

download_path=$1
dir_path=$2
bam_name=$3

if ! scp -P 1300 autosome@localhost:"$download_path" "$dir_path/${bam_name}"
then
  echo "Failed to download $download_path in $dir_path/${bam_name}"
  exit 1
fi

mv "$dir_path/${bam_name}" /mnt/NAS/home/abramov/atac/



