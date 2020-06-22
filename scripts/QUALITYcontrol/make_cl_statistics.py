import gzip
import os.path
import sys
import pandas as pd

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.helpers import remove_punctuation, make_list_from_vcf_without_filter, make_list_from_vcf
from scripts.HELPERS.paths import create_path_from_GTRD_function
from scripts.HELPERS.paths_for_components import GTRD_slice_path, parameters_path

# interestingSet = {'H1 (embryonic stem cells)',
#                   'BGO3 (embryonic stem cells)',
#                   'HUES64 (embryonic stem cells)',
#                   'CyT49 (embryonic stem cells)',
#                   'H9 (embryonic stem cells)',
#                   'VAL-3 (embryonic stem cells)',
#                   'WA09 (embryonic stem cells)',
#                   'embryonic stem cells',
#                   'human embryonic stem cells, H1 (WA01)'}

# interestingSet = {"K562 (myelogenous leukemia)"}
#
# interestingSet = {"HCT-116 (colon carcinoma)"}
# interestingSet = {remove_punctuation(x) for x in interestingSet}
SNP_statistics_dict = {}
tf_set = set()
vcf_counter = 0
counter = 0
cl_set = set()
counted_ctrls = set()
ctrl_snp_counter = 0
num = 0
with open(GTRD_slice_path, "r") as ml:
    master_list = ml.readlines()
for line in master_list:
    if line[0] == "#":
        continue
    num +=1
    line = line.split("\t")
    if num % 50 == 0:
        print(num)
    # if remove_punctuation(line[4]) not in interestingSet:
    #     continue
    vcf_path = create_path_from_GTRD_function(line, for_what="vcf", ctrl=True)
    if os.path.isfile(vcf_path) and vcf_path not in counted_ctrls:
        counted_ctrls.add(vcf_path)
        vcf_counter += 1
        with gzip.open(vcf_path, "rt") as vcf_buffer:
            list_of_snps = make_list_from_vcf(vcf_buffer, filter_no_rs=True)
            ctrl_snp_counter += len(list_of_snps)
            for chr, pos, rs_id, ref, alt, ref_counts, alt_counts in list_of_snps:
                try:
                    SNP_statistics_dict[(ref_counts, alt_counts)] += 1
                except KeyError:
                    SNP_statistics_dict[(ref_counts, alt_counts)] = 1

    vcf_path = create_path_from_GTRD_function(line, for_what="vcf")
    if not os.path.isfile(vcf_path):
        continue
    vcf_counter +=1
    tf_set.add(line[1])
    cl_set.add(line[4])
    with gzip.open(vcf_path, "rt") as vcf_buffer:
        list_of_snps = make_list_from_vcf(vcf_buffer, filter_no_rs=True)
        counter += len(list_of_snps)
        for chr, pos, rs_id, ref, alt, ref_counts, alt_counts in list_of_snps:
            try:
                SNP_statistics_dict[(ref_counts, alt_counts)] += 1
            except KeyError:
                SNP_statistics_dict[(ref_counts, alt_counts)] = 1
df = pd.DataFrame({'ref': [], 'alt': [], 'count': []})
for ref, alt in SNP_statistics_dict:
    df = df.append(pd.DataFrame({'ref': [ref], 'alt': [alt], 'count': [SNP_statistics_dict[(ref, alt)]]}))
df.to_csv(os.path.join(parameters_path, "all_snps_statistics.tsv"), sep="\t", index=False)

print('Total snp calls {}, different TFs {}, different cell types {}, VCFs {}'.format(counter + ctrl_snp_counter,
      len(tf_set), len(cl_set), vcf_counter))

