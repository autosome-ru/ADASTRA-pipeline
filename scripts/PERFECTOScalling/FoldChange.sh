#!/bin/bash

source ./Config.cfg

GETNAME(){
	local var=$1
	local varpath=${var%/*}
	[ "$varpath" != "$var" ] && local vartmp="${var:${#varpath}}"
		echo "${vartmp%_*_*}"
}

REFERENCE="/home/abramov/REFERENCE/"
PWM="/home/abramov/PERFECTOScalling/pwms"
FA=$REFERENCE/"genome-norm.fasta"
OUT="/home/abramov/RESULTS/TF_FC/FilteredMaxCover/"
for file in /home/abramov/RESULTS/TFs_for_PERFECTOS/*
do
   
	EXPFILE=$( GETNAME "$file" )
	EXPNAME=${EXPFILE%.*}

	echo $PWM"$EXPNAME"/
	if [ -d $PWM/"$EXPNAME"/ ]; then
		# shellcheck disable=SC2154
		$python3 extract_ape_data.py "$file" $FA "${OUT}${EXPNAME}_ape_data.txt"

		if [ $? != 0 ]; then
    			echo "Failed to extract adjacent nucleotides"
    			continue
		fi

		if [ -f "${OUT}${EXPNAME}_ape_data.txt" ]; then
			echo "Make perfectos"
			# shellcheck disable=SC2154
			$Java -cp ape.jar ru.autosome.perfectosape.SNPScan $PWM/"$EXPNAME"/ "${OUT}${EXPNAME}_ape_data.txt" -P 1 -F 1 > "${OUT}${EXPNAME}_ape.txt"

			if [ $? != 0 ]; then
    				echo "Failed perfectos-ape"
    				continue
			fi

			$python3 adjust_table.py $file "${OUT}${EXPNAME}_ape.txt" "${OUT}${EXPNAME}_fc.txt"

			if [ $? != 0 ]; then
				echo "Failed to add fc to the table"
				continue
			fi

			rm "${OUT}${EXPNAME}_ape_data.txt"
			rm "${OUT}${EXPNAME}_ape.txt"
		else
			echo "NO ASB found for ${EXPNAME}"
		fi
	else
		echo "No PWM found for $EXPNAME"
	fi
done
