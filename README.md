# ADASTra pipeline release Soos 2020-07-14
A pipeline for processing ChIP-seq read allignments in _bam_ format to find allele-specific TF binding (ASB) events.
It consists of 5 main parts:
### A. SNP calling
This part uses [GATK](https://github.com/broadinstitute/gatk/releases) and [PICARD](https://broadinstitute.github.io/picard/) tools for variant calling.
The result is vcf file of SNV calls in [GATK vcf format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format#:~:text=VCF%2C%20or%20Variant%20Call%20Format,indel%2C%20and%20structural%20variation%20calls.&text=This%20document%20describes%20the%20key,output%20by%20the%20GATK%20tools.).
### B. Peak annotation and filtering
Homozygous SNVs, SNVs with less than 5 reads on each allele and not present in dbSNP common colection are filtered out from the vcf files obtained on the previous step. The resulting variants are annotated with ChIP-seq peaks from 4 different peak callers (if available in _bed_ format).
### C. BAD calling
Background Allelic Dosag (BAD) estimation and full-genome BAD maps construction.
See [BABACHI](https://github.com/autosome-ru/BABACHI).
### D. Negative Binomial Mixture fit
Fitting read count distributions separately for reference and alternative alleles and each BAD with Negative Binomial Mixtures.
### E. Statistical evaluation of ASB
Performing one-tailed tests and aggregating the resulting P-values on TF and cell type level using Mudholkar-George method, FDR-correcting the resulting P-values. Evaluating ASB Effect Size.

## Execution and installation
1. Clone this repository to your machine or server
```
git clone https://github.com/autosome-ru/ADASTRA-pipeline/
```
2. Fill the paths to the required files (listed below) in Configs/config.cfg file.
3. From /scripts/ execute pipline_start.sh <stage>
stage is a flag, corresponding to a part of pipeline you wish to start with (listed in order):
- --create_reference create normalized genome and index
- --snp_call GATK snp calling
- --peak_call peak annotation and filtering
- --bad_groups download and annotate BAD groups by GEO GSE or ENCODE id
- --bad_call BAD estimation
- --nb_fit fit negative binomial distributions
- --p_value_count evaluate statistical significance
- --aggregate_p_values perform cell-type and TF-level aggregation of p-values
## Required software
### General
1. Java SE 8
2. Python >= 3.6
3. GATK >= 4.0.12.0
4. PICARD
5. Parallel
### Python packages
numpy>=1.19.0 <br>
pandas>=1.1.0 <br>
scipy>=1.5.1 <br>
statsmodels>=0.11.1 <br>

## Required files
To run the pipeline successfully one must ensure that following files and directories are created in advance and fill absolute path for each file in the scripts/Configs/config.cfg file.
### Directories
- alignments_path = "/home/user/Alignments/"
The directory with .bam files of experiment and control alignments. Should contain 2 subfolders: CTRL/ and EXP/, EXP/ dir should additionaly contain a subfolder for every TF (e.g. EXP/CTCF_HUMAN). Each TF directory and the CTRL/ directory should contain a subdir for every experiment with this TF (e.g. EXP/CTCF_HUMAN/EXP002547/), every one of which should contain exactly one .bam file.
- badmaps_path = "/home/user/BADMAPS/"
A directory to save badmaps in.
- parallel_parameters_path = "/home/user/parallel/"
A directory to save temporary parameters files
- results_path = "/home/user/DATA/"
A directory to save final ASB calls into.
- reference_path = "/home/user/REFERENCE/"
A directory to put normalized genome into. 
- intervals_path = "/home/user/interval/"
A directory with peak calling data. Should contain a subdir for every caller (e.g. MACS), in each of which should be zipped bed-like files with peak calls (names are arbitrary, ending with .interval.zip). However, peaks from different callers, but for the same experiment must have the same name.

### Files
- master_list_path = "/home/user/PARAMETERS/Master-lines.tsv" <br>
A .tsv file with the following columns, each row corresponding to a single experiment:
id - Unique experiment identifier. Must correspond to the experiment folder in alignments_path/EXP/<TF_NAME>/ or alignments_path/CTRL/ if is input.
tf_uniprot_id - TF uniprot name, e.g. CTCF_HUMAN (or arbitrary TF identifier). Must correspond to the name of the corresponding TF directory in alignment_path/EXP/ <br>
species - unused <br>
cell_id - unused <br>
cell_title - Name or identifier of cell type. Used in BADmaps groupping. <br>
READS - unused <br>
ALIGNS - name of corresponding .bam file without extention ('.bam'). <br>
PEAKS - name of corresponding peak call files (without .interval.zip) or 'None' <br>
GEO_GSE - GSE of the study or 'None' <br>
ENCODE - encode id of the experiment or 'None' <br>
control_id - same for corresponding control experiment or nothing if the experiment is control itself (for all below) <br>
control_cell_id <br>
control_cell_title <br>
control_READS <br>
control_ALIGNS <br>
control_GEO_GSE <br>
control_ENCODE <br>
- genome_path = "/home/user/REFERENCE/genome.fasta"
Path to the reference genome file.
- dbsnp_vcf_path = "/home/user/REFERENCE/dbsnp_common.vcf.gz"
Path to dbsnp common collection (gzipped)
- repeats_path = "/home/user/repeats"
Path to repeat annotation .bed file



