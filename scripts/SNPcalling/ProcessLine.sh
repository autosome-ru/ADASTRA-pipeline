#!/bin/bash

AlignmentsPath="/home/abramov/Alignments/"
ScriptsPath="/home/abramov/ASB-Project/scripts/"
SNPcallingScriptsPath=${ScriptsPath}"SNPcalling/"


LINE=$1
IFS=$'\t'
read -ra ADDR <<< "$LINE"

ExpName=${ADDR[0]}
TF=${ADDR[1]}
ReadGroups=${ADDR[5]}
AlignName=${ADDR[6]}
AlignmentDownloadPath=${ADDR[7]}

if [ "$TF" != "None" ]; then
  if [ -f ${AlignmentsPath}"EXP/$TF/$ExpName/$AlignName.vcf" ];then
    rm ${AlignmentsPath}"EXP/$TF/$ExpName/$AlignName.vcf"
  fi
  AlignmentFullPath=${AlignmentsPath}"EXP/$TF/$ExpName/$AlignName.bam"
  OutPath=${AlignmentsPath}"EXP/$TF/$ExpName/"
else
  if [ -f ${AlignmentsPath}"CTRL/$ExpName/$AlignName.vcf.gz" ];then
    rm ${AlignmentsPath}"CTRL/$ExpName/$AlignName.vcf.gz"
  fi
  AlignmentFullPath=${AlignmentsPath}"CTRL/$ExpName/$AlignName.bam"
  OutPath=${AlignmentsPath}"CTRL/$ExpName/"
fi

bash ${SNPcallingScriptsPath}DownloadBam.sh "$AlignmentDownloadPath" "$AlignmentFullPath"

bash ${SNPcallingScriptsPath}AddReadGroups.sh "$AlignmentFullPath" "$ReadGroups"

echo "Doing SNPcalling for $TF $ExpName"
bash ${SNPcallingScriptsPath}SNPcalling.sh -Exp "$AlignmentFullPath" \
	-Out "$OutPath"

if [ $? != 0 ]; then
  echo "Failed SNPcalling $ExpName"
  exit 1
fi

rm "$AlignmentFullPath"
rm "$AlignmentFullPath.bai"

gzip "$OutPath$AlignName.vcf"

if [ $? != 0 ]; then
	echo "Failed gzip vcf $ExpName"
	exit 1
fi


