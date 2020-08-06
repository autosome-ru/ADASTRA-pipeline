#!/bin/bash

source ../Config/CONFIG.cfg
source ../HELPERS/paths_for_components.py

PEAKannotationScriptsPath=${scripts_path}"/PEAKannotation/"

LINE=$1
IFS=$'\t'
read -ra ADDR <<< "$LINE"
	ExpName=${ADDR[0]}
	TF=${ADDR[1]}
	AlignName=${ADDR[6]}
	PeaksName=${ADDR[7]}
if [ "$ExpName" == "#*" ]; then
  exit 1
fi
VCFPath="${alignments_path}/EXP/$TF/$ExpName/$AlignName.vcf.gz"
echo "$VCFPath"

echo "Making $ExpName"
echo "Checking exp VCF"
if ! [ -f "$VCFPath" ]; then
  echo "There is no VCF for exp $ExpName ($TF)"
  exit 1
fi

if  [ -f "${alignments_path}/EXP/$TF/$ExpName/${AlignName}_table_annotated.txt" ]; then
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


if ! bash ${PEAKannotationScriptsPath}/MakeAnnotatedTable.sh -Out $alignments_path/EXP/"$TF/$ExpName" \
		-Rep "$repeats_path" \
		$PeakM $PEAKM $PeakS $PEAKS $PeakG $PEAKG $PeakC $PEAKC\
		-VCF "$VCFPath"
then
  echo "Failed to make tables"

fi

