# ADASTra pipeline release Soos 2020-07-14
A pipeline for processing ChIP-seq read allignments in _bam_ format to find allele-specific TF binding (ASB) events.
It consists of 5 main parts:
### 1. SNP calling
This part uses [GATK](https://github.com/broadinstitute/gatk/releases) and [PICARD](https://broadinstitute.github.io/picard/) tools for variant calling.
The result is vcf file of SNV calls in [GATK vcf format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format#:~:text=VCF%2C%20or%20Variant%20Call%20Format,indel%2C%20and%20structural%20variation%20calls.&text=This%20document%20describes%20the%20key,output%20by%20the%20GATK%20tools.).
### 2. Peak annotation and filtering
Homozygous SNVs, SNVs with less than 5 reads on each allele and not present in dbSNP common colection are filtered out from the vcf files obtained on the previous step. The resulting variants are annotated with ChIP-seq peaks from 4 different peak callers (if available in _bed_ format).
### 3. BAD calling
Background Allelic Dosag (BAD) estimation and full-genome BAD maps construction.
See [BABACHI](https://github.com/autosome-ru/BABACHI).
### 4. Negative Binomial Mixture fit
Fitting read count distributions separately for reference and alternative alleles and each BAD with Negative Binomial Mixtures.
### 5. Statistical evaluation of ASB
Performing one-tailed tests and aggregating the resulting P-values on TF and cell type level using Mudholkar-George method, FDR-correcting the resulting P-values. Evaluating ASB Effect Size.

## Required software
### General
1. Java SE 8
2. Python >= 3.6
3. GATK >= 4.0.12.0
4. PICARD
### Python packages
numpy>=1.19.0 <br>
pandas>=1.1.0 <br>
scipy>=1.5.1 <br>
statsmodels>=0.11.1 <br>

## Reqired files
### Directories
alignments_path="/home/user/Alignments/" <br>
parameters_path="/home/user/PARAMETERS/" <br>
ploidy_path="/home/user/BAD/" <br>
results_path="/home/user/DATA/" <br>
reference_path="/home/user/REFERENCE/" <br>
### Files
data_slice_path="/home/user/PARAMETERS/Master-lines.tsv" <br>
GENOME="/home/user/REFERENCE/genome.fasta" <br>
DBSNP_VCF="/home/user/REFERENCE/dbsnp_common.vcf.gz" <br>


