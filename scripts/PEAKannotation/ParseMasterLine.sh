#!/bin/bash

source ../../CONFIG.cfg
source ../HELPERS/paths_for_components.py

PEAKannotationScriptsPath=${scripts_path}"/PEAKannotation/"

IFS=$'\t'
read -ra ADDR <<< "$1"

base_path=${ADDR[0]}
peaks_name=${ADDR[1]}

vcf_path="${base_path}.vcf.gz"
echo "$vcf_path"

echo "Making ${base_path}"
echo "Checking exp VCF"
if ! [ -f "$vcf_path" ]; then
  echo "There is no VCF for exp $vcf_path"
  exit 1
fi

if  [ -f "${base_path}_table_annotated.tsv" ]; then
	echo "Remaking annotation for $vcf_path"
else
	echo "Making annotation for $vcf_path first time"
fi

if [ -f "$intervals_path/macs/${peaks_name}.interval.zip" ];then
	PeakM="-macs"
	PEAKM="$intervals_path/macs/${peaks_name}.interval.zip"
else
  PeakM=""
  PEAKM=""
fi

if [ -f "$intervals_path/gem/${peaks_name}.interval.zip" ];then
  PeakG="-gem"
  PEAKG="$intervals_path/gem/${peaks_name}.interval.zip"
else
  PeakG=""
  PEAKG=""
fi

if [ -f "$intervals_path/cpics/${peaks_name}.interval.zip" ];then
  PeakC="-cpics"
  PEAKC="$intervals_path/cpics/${peaks_name}.interval.zip"
else
  PeakC=""
  PEAKC=""
fi

if [ -f "$intervals_path/sissrs/${peaks_name}.interval.zip" ];then
  PeakS="-sissrs"
  PEAKS="$intervals_path/sissrs/${peaks_name}.interval.zip"
else
  PeakS=""
  PEAKS=""
fi


if ! bash ${PEAKannotationScriptsPath}/MakeAnnotatedTable.sh -Out ${base_path} \
		-Rep "$repeats_path" \
		$PeakM $PEAKM $PeakS $PEAKS $PeakG $PEAKG $PeakC $PEAKC\
		-VCF "$vcf_path"
then
  echo "Failed to make annotation table for ${vcf_path}"

fi

