#!/bin/bash

source ../HELPERS/Config.cfg

GETNAME(){
	local var=$1
	local varpath=${var%/*}
	[ "$varpath" != "$var" ] && local vartmp="${var:${#varpath}}"
		echo "${vartmp%_*_*}"
}

ReferencePath="/home/abramov/ReferencePath/"
PWMs_path="/home/abramov/PERFECTOScalling/pwms"
FA=$ReferencePath/"genome-norm.fasta"
OutPath="/home/abramov/RESULTS/TF_FC/FilteredMaxCover/"
ThresholdsPath="/home/abramov/ThresholdsPath"
ResultsPath="/home/abramov/RESULTS/"
for file in "${ResultsPath}TFs_for_PERFECTOS/"*
do
	ExpFile=$( GETNAME "$file" )
	ExpName=${ExpFile%.*}

	echo $PWMs_path"$ExpName"/
	if [ -d $PWMs_path/"$ExpName"/ ]; then
		# shellcheck disable=SC2154


		if ! $python3 extract_ape_data.py "$file" $FA "${OutPath}${ExpName}_ape_data.txt"
		then
    			echo "Failed to extract adjacent nucleotides"
    			continue
		fi

		if [ -f "${OutPath}${ExpName}_ape_data.txt" ]; then
			echo "Make perfectos"
			# shellcheck disable=SC2154

			if ! $Java -cp ape.jar ru.autosome.perfectosape.SNPScan $PWMs_path/"$ExpName/" \
			                        "${OutPath}${ExpName}_ape_data.txt" --precalc "$ThresholdsPath" \
			                        -P 1 -F 1 > ${OutPath}
			then
    				echo "Failed perfectos-ape"
    				continue
			fi



			if ! $python3 adjust_table.py $file "${OutPath}${ExpName}_ape.txt" "${OutPath}${ExpName}_fc.txt";
			then
				echo "Failed to add fc to the table"
				continue
			fi

			rm "${OutPath}${ExpName}_ape_data.txt"
			rm "${OutPath}${ExpName}_ape.txt"
		else
			echo "NO ASB found for ${ExpName}"
		fi
	else
		echo "No PWMs_path found for $ExpName"
	fi
done
