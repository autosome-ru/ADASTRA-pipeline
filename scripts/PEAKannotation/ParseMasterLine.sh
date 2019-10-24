#!/bin/bash

AlignmentsPath="/home/abramov/Alignments/"
IntervalsPath="/home/abramov/intervals/"
RepPath="/home/abramov/Repeats/repeats_regions.tsv"
LINE=$1
IFS=$'\t'
read -ra ADDR <<< "$LINE"
	ExpName=${ADDR[0]}
	TF=${ADDR[1]}
	AlignName=${ADDR[6]}
	PeaksName=${ADDR[7]}

echo "${AlignmentsPath}EXP/$TF/$ExpName/$AlignName.vcf.gz"

echo "Making $ExpName"
echo "Checking exp VCF"
if ! [ -f "${AlignmentsPath}EXP/$TF/$ExpName/$AlignName.vcf.gz" ]; then
  echo "There is no VCF for exp $ExpName ($TF)"
  exit
fi

if  [ -f "${AlignmentsPath}EXP/$TF/$ExpName/${AlignName}_table_annotated.txt" ]; then
	echo "Remaking $ExpName"
else
	echo "Making $ExpName first time"
fi

if [ -f "$IntervalsPath/macs/${PeaksName}.interval.zip" ];then
	PeakM="-macs"
	PEAKM="$IntervalsPath/intervals/macs/${PeaksName}.interval.zip"
else
  PeakM=""
  PEAKM=""
fi

if [ -f "$IntervalsPath/gem/${PeaksName}.interval.zip" ];then
  PeakG="-gem"
  PEAKG="$IntervalsPath/gem/${PeaksName}.interval.zip"
else
  PeakG=""
  PEAKG=""
fi

if [ -f "$IntervalsPath/cpics/${PeaksName}.interval.zip" ];then
  PeakC="-cpics"
  PEAKC="$IntervalsPath/cpics/${PeaksName}.interval.zip"
else
  PeakC=""
  PEAKC=""
fi

if [ -f "$IntervalsPath/sissrs/${PeaksName}.interval.zip" ];then
  PeakS="-sissrs"
  PEAKS="$IntervalsPath/sissrs/${PeaksName}.interval.zip"
else
  PeakS=""
  PEAKS=""
fi

bash MakeAnnotatedTable.sh -Out $AlignmentsPath/ExpName/"$TF/$ExpName" \
		-Rep "$RepPath" \
		$PeakM $PEAKM $PeakS $PEAKS $PeakG $PEAKG $PeakC $PEAKC\
		-VCF $AlignmentsPath/ExpName/"$TF/$ExpName/$AlignName.vcf.gz"
if [ $? != 0 ]; then
  echo "Failed to make tables"

fi

