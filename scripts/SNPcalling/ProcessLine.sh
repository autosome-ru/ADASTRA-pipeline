#!/bin/bash

source ../HELPERS/paths_for_components.py
source ../Configs/CONFIG.cfg

snp_calling_scripts_path=${scripts_path}"/SNPcalling/"
helpers_scripts_path=${scripts_path}"/HELPERS/"

LINE=$2
to_download=$1
IFS=$'\t'
read -ra ADDR <<< "$LINE"

ExpName=${ADDR[0]}
ReadGroups=${ADDR[5]}
AlignName=${ADDR[6]}
AlignmentDownloadPath=${ADDR[7]}

if [ "$to_download" == "-d" ]; then
  if [ "$AlignmentDownloadPath" = "None" ];then
      echo "There is no Path for exp $ExpName"
      exit 1
  fi
fi

echo "Making dirs"
if ! [ -d ${alignments_path}"$ExpName" ]; then
  if ! mkdir ${alignments_path}"$ExpName"
  then
    echo "Failed to make dir $ExpName"
    exit 1
  fi
fi

OutPath=${alignments_path}"$ExpName/"
AlignmentFullPath=${alignments_path}"/$ExpName/$AlignName.bam"


echo "Downloading $ExpName"
if [ "$to_download" == "-d" ]; then
  if ! bash ${helpers_scripts_path}DownloadFile.sh "$AlignmentDownloadPath" "$AlignmentFullPath"
  then
    echo "Download failed for $ExpName"
    exit 1
  fi
fi

echo "Adding ReadGroups for $ExpName"
if ! bash ${snp_calling_scripts_path}AddReadGroups.sh "$AlignmentFullPath" "$ReadGroups"
then
  echo "Failed AddReadGroups $ExpName"
  exit 1
fi

echo "Doing SNPcalling for $TF $ExpName"
if ! bash ${snp_calling_scripts_path}SNPcalling.sh -Exp "$AlignmentFullPath" -Out "$OutPath"
then
  echo "Failed SNPcalling $ExpName"
  exit 1
fi

rm "$AlignmentFullPath"
rm "$AlignmentFullPath.bai"

if [ -f ${OutPath}"$AlignName.vcf.idx" ];then
  rm ${OutPath}"$AlignName.vcf.idx"
fi

if [ -f ${OutPath}"$AlignName.vcf.gz" ];then
  rm ${OutPath}"$AlignName.vcf.gz"
fi

if ! gzip "$OutPath$AlignName.vcf"
then
	echo "Failed gzip vcf $ExpName"
	exit 1
fi


