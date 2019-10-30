#!/bin/bash
source ../HELPERS/Config.cfg

VCF="/home/abramov/REFERENCE/00-common_all.vcf.gz"
REFERENCE="/home/abramov/REFERENCE/"
FA=${REFERENCE}"genome-norm.fasta"

while [ "$(echo "$1" | cut -c1)" = "-" ]
do
    case "$1" in
    -Out) OutPath=$2
        shift 2;;
	  -Exp) BAM=$2
		    BamPath=${BAM%/*}
		    [ "$BamPath" != "$BAM" ] && TMP="${BAM:${#BamPath}}"
        BamName=${TMP%.*}
        shift 2;;
    *)
        echo "There is no option $1"
	      break;;
	  esac
done

bam_size=0

if [ ! -f "$VCF.tbi" ]; then
	echo "Index file for VCF not found, indexing.."
	# shellcheck disable=SC2154
	$Java $JavaParameters  -jar "$GATK" \
		IndexFeatureFile \
       		"-F $VCF"
fi

if [  -f "$BamPath/$BamName.bam.bai" ]; then
	rm "$BamPath/$BamName.bam.bai"
fi	
echo "Index file for $BamPath/$BamName"

if ! samtools index "$BamPath/$BamName.bam"
then
        echo "Failed to index bam"
        exit 1
fi

echo "Cutting bam.."
if ! samtools view -b "$BamPath/$BamName.bam" \
	chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 \
	chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > "$OutPath/${BamName}_chop.bam"
then
    echo "Failed to cut bam"
    exit 1
fi

if ! $Java $JavaParameters -jar "$PICARD" \
	AddOrReplaceReadGroups \
	I="$OutPath/${BamName}_chop.bam" \
	O="$OutPath/${BamName}_formated.bam" \
	RGID=1 \
	RGLB=lib1 \
	RGPL=seq1 \
	RGPU=unit1 \
	RGSM=20
then
    echo "Failed to picard.AddOrReplaceReadGroups"
    exit 1
fi

bam_size=$((bam_size + $(wc -c <"$OutPath${BamName}_chop.bam")))
rm "$OutPath${BamName}_chop.bam"

if ! $Java $JavaParameters -jar "$PICARD" \
	MarkDuplicates \
	I="$OutPath/${BamName}_formated.bam" \
	O="$OutPath/${BamName}_ready.bam" \
	REMOVE_DUPLICATES=true \
	M="$OutPath/${BamName}_metrics.txt"
then
    echo "Failed to picard.MarkDplicates"
    exit 1
fi

bam_size=$((bam_size + $(wc -c <"$OutPath${BamName}_formated.bam")))
rm "$OutPath${BamName}_formated.bam"
rm "$OutPath/${BamName}_metrics.txt"

if ! $Java $JavaParameters -jar "$GATK" \
	BaseRecalibrator \
	-R "$FA" \
	-I "$OutPath/${BamName}_ready.bam" \
	-known-sites "$VCF" \
	-O "$OutPath/${BamName}.table"
then
    echo "Failed to make base recalibration"
    exit 1
fi

if ! $Java $JavaParameters -jar "$GATK" \
	ApplyBQSR \
	-R "$FA" \
	-I "$OutPath/${BamName}_ready.bam" \
	--bqsr-recal-file "$OutPath/${BamName}.table" \
	-O "$OutPath/${BamName}_final.bam"
then
    echo "Failed to apply BSQR"
    exit 1
fi

bam_size=$((bam_size + $(wc -c <"$OutPath${BamName}_ready.bam")))
if [ -f "$OutPath/${BamName}_ready.bam.bai" ]; then
    rm "$OutPath/${BamName}_ready.bam.bai"
fi
rm "$OutPath/${BamName}.table"
rm "$OutPath${BamName}_ready.bam"

if ! $Java $JavaParameters -jar "$GATK" \
	HaplotypeCaller \
	-R "$FA" \
	-I "$OutPath/${BamName}_final.bam" \
	--dbsnp "$VCF" \
	-O "$OutPath/${BamName}.vcf"
then
    echo "Failed gatk.HaplotypeCaller"
    exit 1
fi

bam_size=$((bam_size + $(wc -c <"$OutPath${BamName}_final.bam")))

rm "$OutPath${BamName}_final.bam"
rm "$OutPath${BamName}_final.bai"

echo "Total intermediate .bam size: $bam_size"
