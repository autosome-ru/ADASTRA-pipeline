#!/bin/bash

echo ../../CONFIG.cfg
source ../../CONFIG.cfg


while [ "$(echo "$1" | cut -c1)" = "-" ]
do
  case "$1" in
    RefFolder) OUT=$2
      shift 2;;

	  RefGenome) REF=$2
      shift 2;;
    *)
      echo "There is no option $1"
		  break;;
	esac
done

if ! $Java $JavaParameters -jar "$PICARD" \
	NormalizeFasta \
	I="$REF" \
	O="$OUT/genome-norm.fasta"
then
    echo "Failed to normalize fasta"
    exit 1
fi

if ! samtools faidx "$OUT/genome-norm.fasta"
then
    echo "Failed to index fasta"
    exit 1
fi

if ! $Java $JavaParameters -jar "$PICARD" \
	CreateSequenceDictionary \
	R="$OUT/genome-norm.fasta"\
	O="$OUT/genome-norm.dict"
then
    echo "Failed to create sequence dictionary"
    exit 1
fi