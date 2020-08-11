#!/bin/bash

samtools view -H "$1" > "$1.header.txt"
IFS=","
read -ra A <<< "$2"
if [ ${#A[@]} -le 1 ]; then
  exit 0
fi
for head in "${A[@]}";do
	echo "$head"
	echo -e "@RG\tID:$head\tPL:ILLUMINA\tLB:library\tSM:sample\tPI:400" >> "$1.header.txt"
done
samtools reheader "$1.header.txt" "$1" > "$1.new"
rm "$1"
mv "$1.new" "$1"
rm "$1.header.txt"
echo "$1"