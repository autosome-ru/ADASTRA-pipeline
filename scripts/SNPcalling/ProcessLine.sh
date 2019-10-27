#!/bin/bash

AlignmentsPath="/home/abramov/Alignments/"
ScriptsPath="/home/abramov/ASB-Project/scripts/"
SNPcallingScriptsPath=${ScriptsPath}"SNPcalling/"
HelpersScriptsPath=${ScriptsPath}"HELPERS/"

LINE=$2
to_download=$1
IFS=$'\t'
read -ra ADDR <<< "$LINE"

ExpName=${ADDR[0]}
TF=${ADDR[1]}
ReadGroups=${ADDR[5]}
AlignName=${ADDR[6]}
AlignmentDownloadPath=${ADDR[7]}

if [ "$AlignmentDownloadPath" = "None" ];then
    echo "There is no Path for exp $ExpName"
    exit 1
fi

if [ "$TF" != "None" ]; then
  if ! [ -d ${AlignmentsPath}"EXP/$TF" ]; then
    if ! mkdir ${AlignmentsPath}"EXP/$TF"
    then
      echo "Failed to make dir $TF"
      exit 1
    fi
  fi

  if ! [ -d ${AlignmentsPath}"EXP/$TF/$ExpName" ]; then
    if ! mkdir ${AlignmentsPath}"EXP/$TF/$ExpName"
    then
      echo "Failed to make dir $ExpName"
      exit 1
    fi
  else
    echo "Directory for $ExpName already exists"
  fi

  if [ -f ${AlignmentsPath}"EXP/$TF/$ExpName/$AlignName.vcf" ];then
    rm ${AlignmentsPath}"EXP/$TF/$ExpName/$AlignName.vcf"
  fi
  AlignmentFullPath=${AlignmentsPath}"EXP/$TF/$ExpName/$AlignName.bam"
  OutPath=${AlignmentsPath}"EXP/$TF/$ExpName/"
else
  if ! [ -d /home/abramov/Alignments/CTRL/"$EXP" ]; then
    if ! mkdir /home/abramov/Alignments/CTRL/"$EXP"
    then
      echo "Failed to make dir $ExpName"
      exit 1
    fi
  else
    echo "Directory for $ExpName already exists"
  fi

  if [ -f ${AlignmentsPath}"CTRL/$ExpName/$AlignName.vcf.gz" ];then
    rm ${AlignmentsPath}"CTRL/$ExpName/$AlignName.vcf.gz"
  fi
  AlignmentFullPath=${AlignmentsPath}"CTRL/$ExpName/$AlignName.bam"
  OutPath=${AlignmentsPath}"CTRL/$ExpName/"
fi
echo "Downloading $ExpName"
if [ "$to_download" == "-d" ]; then
  if ! bash ${HelpersScriptsPath}DownloadFile.sh "$AlignmentDownloadPath" "$AlignmentFullPath"
  then
    echo "Download failed for $ExpName"
    exit 1
  fi
fi

echo "Adding ReadGroups for $ExpName"
if ! bash ${SNPcallingScriptsPath}AddReadGroups.sh "$AlignmentFullPath" "$ReadGroups"
then
  echo "Failed AddReadGroups $ExpName"
  exit 1
fi

echo "Doing SNPcalling for $TF $ExpName"
if ! bash ${SNPcallingScriptsPath}SNPcalling.sh -Exp "$AlignmentFullPath" -Out "$OutPath"
then
  echo "Failed SNPcalling $ExpName"
  exit 1
fi

rm "$AlignmentFullPath"
rm "$AlignmentFullPath.bai"

if ! gzip "$OutPath$AlignName.vcf"
then
	echo "Failed gzip vcf $ExpName"
	exit 1
fi


