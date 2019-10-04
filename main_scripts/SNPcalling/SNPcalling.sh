#!/bin/bash
source ./Config.cfg
GETNAME(){
	local var=$1
	local varpath=${var%/*}
	[ "$varpath" != "$var" ] && local vartmp="${var:${#varpath}}"
		echo ${vartmp%.*}
}

WG=false
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
	-Exp) EXP=$2
		EXPPATH=${EXP%/*}
		[ "$EXPPATH" != "$EXP" ] && TMP="${EXP:${#EXPPATH}}"
                EXPNAME=${TMP%.*}
		echo "THIS is $EXPNAME and this is $EXPPATH"
        	shift 2;;
	-VCF) VCF=$2
              	shift 2;;
	-WG) WG=true
		shift 1;;
    *)
        echo "There is no option $1"
	break;;
	esac
done

if [ ! -f "$VCF.tbi" ]; then	
	echo "Index file for VCF not found, indexing.."
	$Java $JavaParameters  -jar $GATK \
		IndexFeatureFile \
       		-F $VCF
fi


FA=$REFERENCE/"genome-norm.fasta"
FD=$REFERENCE/"genome-norm.dict"

bash pre-process.sh $EXPNAME \
	$EXPPATH \
	$OUT \
	$VCF \
	$FA \
	$FD \
	$WG

if [ $? != 0 ]; then
    echo "Failed to pre-process exp"
    exit 1
fi

bam_size=0

bam_size=$(($bam_size+$(wc -c <"$OUT${EXPNAME}_formated.bam")))
bam_size=$(($bam_size+$(wc -c <"$OUT${EXPNAME}_ready.bam")))
bam_size=$(($bam_size+$(wc -c <"$OUT${EXPNAME}_chop.bam")))
bam_size=$(($bam_size+$(wc -c <"$OUT${EXPNAME}_final.bam")))

rm "$OUT${EXPNAME}_final.bam"
rm "$OUT${EXPNAME}_final.bai"
rm "$OUT${EXPNAME}_chop.bam"
rm "$OUT${EXPNAME}_ready.bam"
rm "$OUT${EXPNAME}_formated.bam"

echo "Total intermediate .bam size: $bam_size"

exit 0
