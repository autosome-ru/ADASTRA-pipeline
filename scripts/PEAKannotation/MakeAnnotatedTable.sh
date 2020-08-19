#!/bin/bash
source scripts/HELPERS/soft_configs.cfg
source scripts/HELPERS/paths_for_components.py

PEAKannotationScriptsPath=${scripts_path}/PEAKannotation/

GETNAME(){
	local var=$1
	local varpath=${var%/*}
	[ "$varpath" != "$var" ] && local vartmp="${var:${#varpath}}"
		echo "${vartmp%.*}"
}


withmacs=false
withsissrs=false
withcpics=false
withgem=false

while [ "$(echo "$1" | cut -c1)" = "-" ]
do
  case "$1" in
	  -macs) withmacs=true
		  macs=$2
		  shift 2;;

	  -sissrs) withsissrs=true
		  sissrs=$2
		  shift 2;;

	  -cpics) withcpics=true
		  cpics=$2
		  shift 2;;

	  -gem) withgem=true
		  gem=$2
		  shift 2;;

	  -I) base_path=$2
		  shift 2;;

	  *)
		  echo "There is no option $1"
		  break;;

  esac
done

if [ $withgem != false ]; then
	adastra check_pos_peaks --peak "$gem" --out "${base_path}.gem.bed" --type 'gem'
	# shellcheck disable=SC2154
	if ! bedtools sort -i "${base_path}.gem.bed" > "${base_path}.gem.bed.sorted"
	then
		echo "Failed to sort gem peaks"
		exit 1
	fi
	rm "${base_path}.gem.bed"
fi

if [ $withmacs != false ]; then
	adastra check_pos_peaks --peak "$macs" --out "${base_path}.macs.bed" --type 'macs'

	if ! bedtools sort -i "${base_path}.macs.bed" > "${base_path}.macs.bed.sorted"
	then
		echo "Failed to sort macs peaks"
		exit 1
	fi
  rm "${base_path}.macs.bed"
fi

if [ $withsissrs != false ]; then
	adastra check_pos_peaks --peak "$sissrs" --out "${base_path}.sissrs.bed" --type 'sissrs'

	if ! bedtools sort -i "${base_path}.sissrs.bed" > "${base_path}.sissrs.bed.sorted"
	then
		echo "Failed to sort sissrs peaks"
		exit 1
	fi
  rm "${base_path}.sissrs.bed"
fi

if [ $withcpics != false ]; then
	adastra check_pos_peaks --peak "$cpics" --out "${base_path}.cpics.bed" --type 'cpics'

	if ! bedtools sort -i "${base_path}.cpics.bed" > "${base_path}.cpics.bed.sorted"
	then
		echo "Failed to sort cpics peaks"
		exit 1
	fi
  rm "${base_path}.cpics.bed"
fi

adastra annotate_peaks --base "$base_path"


if [ "$withgem" != false ]; then
	rm "${base_path}.gem.bed.sorted"
fi

if [ "$withcpics" != false ]; then
	rm "${base_path}.cpics.bed.sorted"
fi

if [ "$withmacs" != false ]; then
  rm "${base_path}.macs.bed.sorted"
fi

if [ "$withsissrs" != false ]; then
  rm "${base_path}.sissrs.bed.sorted"
fi
