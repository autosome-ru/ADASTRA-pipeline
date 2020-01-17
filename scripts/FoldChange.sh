#!/bin/bash

source HELPERS/Config.cfg
source HELPERS/paths_for_components.py

GETNAME(){
	local var=$1
	local varpath=${var%/*}
	[ "$varpath" != "$var" ] && local vartmp="${var:${#varpath}}"
		echo "${vartmp%_*_*}"
}

for file in "${results_path}TF_P-values/"*
do
  echo "${results_path}TF_P-values/"
	ExpFile=$( GETNAME "$file" )
	ExpName=${ExpFile%.*}

	echo $PWMs_path"$ExpName"/
	if [ -d $PWMs_path/"$ExpName"/ ]; then
		# shellcheck disable=SC2154

		if ! $python3 "${scripts_path}PERFECTOScalling/"extract_ape_data.py "$file" $FA "${perfectos_path}${ExpName}_ape_data.txt"
		then
    			echo "Failed to extract adjacent nucleotides"
    			continue
		fi

		if [ -f "${perfectos_path}${ExpName}_ape_data.txt" ]; then
			echo "Make perfectos"
			# shellcheck disable=SC2154

			if ! $Java -cp ape.jar ru.autosome.perfectosape.SNPScan $PWMs_path/"$ExpName/" \
			                        "${perfectos_path}${ExpName}_ape_data.txt" \
			                        -P 1 -F 1 > "${perfectos_path}${ExpName}_ape.txt"
			then
    				echo "Failed perfectos-ape"
    				continue
			fi



			if ! $python3 "${scripts_path}PERFECTOScalling/"adjust_table.py $file "${perfectos_path}${ExpName}_ape.txt" "${perfectos_path}${ExpName}_fc.txt";
			then
				echo "Failed to add fc to the table"
				continue
			fi

			rm "${perfectos_path}${ExpName}_ape_data.txt"
			rm "${perfectos_path}${ExpName}_ape.txt"
		else
			echo "NO ASB found for ${ExpName}"
		fi
	else
		echo "No PWMs_path found for $ExpName"
	fi
done
