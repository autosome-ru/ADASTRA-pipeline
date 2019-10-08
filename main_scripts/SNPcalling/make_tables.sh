#!/bin/bash
source ./Config.cfg

GETNAME(){
	local var=$1
	local varpath=${var%/*}
	[ "$varpath" != "$var" ] && local vartmp="${var:${#varpath}}"
		echo ${vartmp%.*}
}

withmacs=false
withsissrs=false
withcpics=false
withgem=false
macs=-1
sissrs=-1
cpics=-1
gem=-1

while [ "`echo $1 | cut -c1`" = "-" ]
do
    case "$1" in
	-Out) OUT=$2
		shift 2;;

	-Ref) REFERENCE=$2
		shift 2;;

	-macs) withmacs=true
		macs=$2
		NAMEM=$(GETNAME $macs)
		shift 2;;

	-sissrs) withsissrs=true
		sissrs=$2
		NAMES=$( GETNAME $sissrs)
		shift 2;;

	-cpics) withcpics=true
		cpics=$2
		NAMEC=$( GETNAME $cpics)
		shift 2;;

	-gem) withgem=true
		gem=$2
		NAMEG=$( GETNAME $gem)
		shift 2;;

	-VCFexp) VCFexp=$2
		tmp=$( GETNAME $VCFexp )
		EXPNAME=${tmp%.*}
		shift 2;;

	*)
		echo "There is no option $1"
		break;;

    esac
done

FA=$REFERENCE/"genome-norm.fasta"
FD=$REFERENCE/"genome-norm.dict"


$python3 Make_tables.py "$OUT${EXPNAME}.vcf.gz" "$OUT${EXPNAME}_table.txt"

if [ $withgem != false ]; then
	$python3 CheckPositive.py $gem
	$Bedtools sort -i $gem > "$gem.sorted"
	if [ $? != 0 ]; then
		echo "Failed to sort gem peaks"
		exit 1
	fi
fi

if [ $withmacs != false ]; then
	$python3 CheckPositive.py $macs
	$Bedtools sort -i $macs > "$macs.sorted"
	if [ $? != 0 ]; then
		echo "Failed to sort macs peaks"
		exit 1
	fi
fi

if [ $withsissrs != false ]; then
	$python3 CheckPositive.py $sissrs
	$Bedtools sort -i $sissrs > "$sissrs.sorted"
	if [ $? != 0 ]; then
		echo "Failed to sort sissrs peaks"
		exit 1
	fi
fi

if [ $withcpics != false ]; then
	$python3 CheckPositive.py $cpics
	$Bedtools sort -i $cpics > "$cpics.sorted"
	if [ $? != 0 ]; then
		echo "Failed to sort cpics peaks"
		exit 1
	fi
fi

$python3 Annotate.py "$OUT${EXPNAME}_table.txt" "$macs.sorted" "$sissrs.sorted" "$cpics.sorted" "$gem.sorted" $withmacs $withsissrs $withcpics $withgem "$OUT${EXPNAME}_table_annotated.txt"

rm "$OUT${EXPNAME}_table.txt"
rm "$OUT${EXPNAME}_table_annotated.txt.m.txt"
rm "$OUT${EXPNAME}_table_annotated.txt.c.txt"
rm "$OUT${EXPNAME}_table_annotated.txt.s.txt"

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
