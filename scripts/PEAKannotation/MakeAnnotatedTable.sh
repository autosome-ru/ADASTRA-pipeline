#!/bin/bash
source ../HELPERS/Config.cfg

GETNAME(){
	local var=$1
	local varpath=${var%/*}
	[ "$varpath" != "$var" ] && local vartmp="${var:${#varpath}}"
		echo "${vartmp%.*}"
}

ScriptsPath="/home/abramov/ASB-Project/scripts/"
PEAKannotationScriptsPath=${ScriptsPath}PEAKannotation/

withmacs=false
withsissrs=false
withcpics=false
withgem=false

while [ "$(echo "$1" | cut -c1)" = "-" ]
do
  case "$1" in

	  -Out) OUT=$2
		  shift 2;;

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

    -Rep) RepFile=$2
      shift 2;;

	  -VCF) VCF=$2
		  tmp=$( GETNAME "$VCF" )
		  EXPNAME=${tmp%.*}
		  shift 2;;

	  *)
		  echo "There is no option $1"
		  break;;

  esac
done

if [ $withgem != false ]; then

	python3 ${PEAKannotationScriptsPath}CheckPositive.py "$gem" "$OUT${EXPNAME}_gem.bed" 'gem'
	# shellcheck disable=SC2154
	bedtools sort -i "$OUT${EXPNAME}_gem.bed" > "$OUT${EXPNAME}_gem.bed.sorted"

	rm "$OUT${EXPNAME}_gem.bed"

	if [ $? != 0 ]; then
		echo "Failed to sort gem peaks"
		exit 1
	fi

fi

if [ $withmacs != false ]; then
	python3 ${PEAKannotationScriptsPath}CheckPositive.py "$macs" "$OUT${EXPNAME}_macs.bed" 'macs'
	bedtools sort -i "$OUT${EXPNAME}_macs.bed" > "$OUT${EXPNAME}_macs.bed.sorted"
  rm "$OUT${EXPNAME}_macs.bed"

	if [ $? != 0 ]; then
		echo "Failed to sort macs peaks"
		exit 1
	fi

fi

if [ $withsissrs != false ]; then
	python3 ${PEAKannotationScriptsPath}CheckPositive.py "$sissrs" "$OUT${EXPNAME}_sissrs.bed" 'sissrs'
	bedtools sort -i "$OUT${EXPNAME}_sissrs.bed" > "$OUT${EXPNAME}_sissrs.bed.sorted"
  rm "$OUT${EXPNAME}_sissrs.bed"

	if [ $? != 0 ]; then
		echo "Failed to sort sissrs peaks"
		exit 1
	fi

fi

if [ $withcpics != false ]; then
	python3 ${PEAKannotationScriptsPath}CheckPositive.py "$cpics" "$OUT${EXPNAME}_cpics.bed" 'cpics'
	bedtools sort -i "$OUT${EXPNAME}_cpics.bed" > "$OUT${EXPNAME}_cpics.bed.sorted"

	rm "$OUT${EXPNAME}_cpics.bed"
	if [ $? != 0 ]; then
		echo "Failed to sort cpics peaks"
		exit 1
	fi

fi

python3 ${PEAKannotationScriptsPath}Annotate.py "$VCF" "$OUT${EXPNAME}_table_annotated.txt" "$RepFile"


if [ "$withgem" != false ]; then
	rm "$OUT${EXPNAME}_gem.bed.sorted"
fi

if [ "$withcpics" != false ]; then
	rm "$OUT${EXPNAME}_cpics.bed.sorted"
fi

if [ "$withmacs" != false ]; then
  rm "$OUT${EXPNAME}_macs.bed.sorted"
fi

if [ "$withsissrs" != false ]; then
  rm "$OUT${EXPNAME}_sissrs.bed.sorted"
fi


exit 0
