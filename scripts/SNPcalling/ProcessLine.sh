#!/bin/bash

source ../HELPERS/paths_for_components.py
source ../Configs/CONFIG.cfg

snp_calling_scripts_path=${scripts_path}"/SNPcalling/"
helpers_scripts_path=${scripts_path}"/HELPERS/"

LINE=$2
to_download=$1
IFS=$'\t'
read -ra ADDR <<< "$LINE"

exp_name=${ADDR[0]}
align_name=${ADDR[1]}
read_groups=${ADDR[2]}
AlignmentDownloadPath=${ADDR[3]}

if [ "$to_download" == "-d" ]; then
  if [ "$AlignmentDownloadPath" = "None" ];then
      echo "There is no Path for exp $exp_name"
      exit 1
  fi
fi

echo "Making dirs"
if ! [ -d ${alignments_path}"$exp_name" ]; then
  if ! mkdir ${alignments_path}"$exp_name"
  then
    echo "Failed to make dir $exp_name"
    exit 1
  fi
fi

OutPath=${alignments_path}"$exp_name/"
AlignmentFullPath=${alignments_path}"/$exp_name/${align_name}.bam"


echo "Downloading $exp_name"
if [ "$to_download" == "-d" ]; then
  if ! bash ${helpers_scripts_path}DownloadFile.sh "$AlignmentDownloadPath" "$AlignmentFullPath"
  then
    echo "Download failed for $exp_name"
    exit 1
  fi
fi

echo "Adding read_groups for $exp_name"
if ! bash ${snp_calling_scripts_path}AddReadGroups.sh "$AlignmentFullPath" "$read_groups"
then
  echo "Failed AddReadGroups $exp_name"
  exit 1
fi

echo "Doing SNPcalling for $TF $exp_name"
if ! bash ${snp_calling_scripts_path}SNPcalling.sh -Exp "$AlignmentFullPath" -Out "$OutPath"
then
  echo "Failed SNPcalling $exp_name"
  exit 1
fi

rm "$AlignmentFullPath"
rm "$AlignmentFullPath.bai"

if [ -f ${OutPath}"$align_name.vcf.idx" ];then
  rm ${OutPath}"$align_name.vcf.idx"
fi

if [ -f ${OutPath}"$align_name.vcf.gz" ];then
  rm ${OutPath}"$align_name.vcf.gz"
fi

if ! gzip "$OutPath$align_name.vcf"
then
	echo "Failed gzip vcf $exp_name"
	exit 1
fi


