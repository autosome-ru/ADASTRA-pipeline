#!/bin/bash
source ./Config.cfg

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
macs=-1
sissrs=-1
cpics=-1
gem=-1

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

	-VCFexp) VCFexp=$2
		tmp=$( GETNAME "$VCFexp" )
		EXPNAME=${tmp%.*}
		shift 2;;

	*)
		echo "There is no option $1"
		break;;

    esac
done


if [ $withgem != false ]; then
	# shellcheck disable=SC2154
	$python3 CheckPositive.py "$gem"
	# shellcheck disable=SC2154
	$Bedtools sort -i "$gem" > "$gem.sorted"

	if [ $? != 0 ]; then
		echo "Failed to sort gem peaks"
		exit 1
	fi

fi

if [ $withmacs != false ]; then
	$python3 CheckPositive.py "$macs"
	$Bedtools sort -i "$macs" > "$macs.sorted"

	if [ $? != 0 ]; then
		echo "Failed to sort macs peaks"
		exit 1
	fi

fi

if [ $withsissrs != false ]; then
	$python3 CheckPositive.py "$sissrs"
	$Bedtools sort -i "$sissrs" > "$sissrs.sorted"

	if [ $? != 0 ]; then
		echo "Failed to sort sissrs peaks"
		exit 1
	fi

fi

if [ $withcpics != false ]; then
	$python3 CheckPositive.py "$cpics"
	$Bedtools sort -i "$cpics" > "$cpics.sorted"

	if [ $? != 0 ]; then
		echo "Failed to sort cpics peaks"
		exit 1
	fi

fi

$python3 Annotate.py "$OUT${EXPNAME}.vcf.gz" "$macs.sorted" "$sissrs.sorted" "$cpics.sorted" "$gem.sorted" \
                      $withmacs $withsissrs $withcpics $withgem "$OUT${EXPNAME}_table_annotated.txt"

if [ "$withgem" != false ]; then
	rm "${gem}.sorted"
fi

if [ "$withcpics" != false ]; then
	rm "${cpics}.sorted"
fi

if [ "$withmacs" != false ]; then
	rm "${macs}.sorted"
fi

if [ "$withsissrs" != false ]; then
	rm "${sissrs}.sorted"
fi


exit 0
