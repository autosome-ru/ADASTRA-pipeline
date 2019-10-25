#!/bin/bash

AlignmentsPath="/home/abramov/Alignments/"

LINE=$1
IFS=$'\t'
read -ra ADDR <<< "$LINE"
if [ "${ADDR[2]}" = "Homo sapiens" ]; then
	ExpName=${ADDR[0]}
	TF=${ADDR[1]}
	AlignName=${ADDR[6]}
fi

if [ "$TF" != "None" ]; then

  if [ -f ${AlignmentsPath}"EXP/$TF/$ExpName/$AlignName.vcf" ];then
    rm ${AlignmentsPath}"EXP/$TF/$ExpName/$AlignName.vcf"
  fi

	echo "Doing SNPcalling for CTRL $ExpName"
	bash SNPcalling.sh -Exp ${AlignmentsPath}"EXP/$TF/$ExpName/$AlignName.bam" \
	-Out ${AlignmentsPath}"EXP/$TF/$ExpName"
	if [ $? != 0 ]; then
    echo "Failed SNPcalling $ExpName"
    exit 1
  fi

  rm ${AlignmentsPath}"EXP/$TF/$ExpName/$AlignName.bam"
  rm ${AlignmentsPath}"EXP/$TF/$ExpName/$AlignName.bam.bai"

	gzip ${AlignmentsPath}"EXP/$TF/$ExpName/$AlignName.vcf"

	if [ $? != 0 ]; then
		echo "Failed SNPcalling $ExpName"
	  exit 1
	fi

else
	if [ -f ${AlignmentsPath}"CTRL/$ExpName/$AlignName.vcf.gz" ];then
    rm ${AlignmentsPath}"CTRL/$ExpName/$AlignName.vcf.gz"
  fi

	echo "Doing SNPcalling for $TF $ExpName"
	bash SNPcalling.sh -Exp ${AlignmentsPath}"CTRL/$ExpName/$AlignName.bam" \
		-Out ${AlignmentsPath}"CTRL/$ExpName"
	if [ $? != 0 ]; then
    echo "Failed SNPcalling $ExpName"
    exit 1
  fi

  rm ${AlignmentsPath}"CTRL/$ExpName/$AlignName.bam"
	rm ${AlignmentsPath}"CTRL/$ExpName/$AlignName.bam.bai"

	gzip ${AlignmentsPath}"CTRL/$ExpName/$AlignName.vcf"

	if [ $? != 0 ]; then
		echo "Failed SNPcalling $ExpName"
	  exit 1
	fi
fi
