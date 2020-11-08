#!/bin/bash

download_path=$1
dir_path=$2
bam_name=$3

if ! scp -T -P 1300 autosome@localhost:"$download_path" "$dir_path"
then
  echo "Failed to download $download_path in $dir_path"
  exit 1
fi

if ! adastra make_bam_list --exp-dir "$dir_path" --download $download_path
then
  echo 'Failed to merge bam'
fi

if ! samtools merge -b "$dir_path/bam_list.txt" "$dir_path"/${bam_name}
then
  echo "Failed to merge vcfs for $download_path"
  exit 1
fi

