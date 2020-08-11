#!/bin/bash

source ../Config/CONFIG.cfg
source ../HELPERS/paths_for_components.py

PEAKannotationScriptsPath=${scripts_path}"/PEAKannotation/"

base_path=$1

vcf_path="${base_path}.vcf.gz"
echo "$vcf_path"

echo "Making ${base_path}"
echo "Checking exp VCF"
if ! [ -f "$vcf_path" ]; then
  echo "There is no VCF for exp $ExpName ($TF)"
  exit 1
fi

if  [ -f "${base_path}_table_annotated.tsv" ]; then
	echo "Remaking $ExpName"
else
	echo "Making $ExpName first time"
fi

if [ -f "$intervals_path/macs/${PeaksName}.interval.zip" ];then
	PeakM="-macs"
	PEAKM="$intervals_path/macs/${PeaksName}.interval.zip"
else
  PeakM=""
  PEAKM=""
fi

if [ -f "$intervals_path/gem/${PeaksName}.interval.zip" ];then
  PeakG="-gem"
  PEAKG="$intervals_path/gem/${PeaksName}.interval.zip"
else
  PeakG=""
  PEAKG=""
fi

if [ -f "$intervals_path/cpics/${PeaksName}.interval.zip" ];then
  PeakC="-cpics"
  PEAKC="$intervals_path/cpics/${PeaksName}.interval.zip"
else
  PeakC=""
  PEAKC=""
fi

if [ -f "$intervals_path/sissrs/${PeaksName}.interval.zip" ];then
  PeakS="-sissrs"
  PEAKS="$intervals_path/sissrs/${PeaksName}.interval.zip"
else
  PeakS=""
  PEAKS=""
fi


if ! bash ${PEAKannotationScriptsPath}/MakeAnnotatedTable.sh -Out ${base_path} \
		-Rep "$repeats_path" \
		$PeakM $PEAKM $PeakS $PEAKS $PeakG $PEAKG $PeakC $PEAKC\
		-VCF "$vcf_path"
then
  echo "Failed to make tables"

fi

