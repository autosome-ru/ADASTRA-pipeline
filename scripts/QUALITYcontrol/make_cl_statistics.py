import os.path
import sys
import pandas as pd

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.helpers import remove_punctuation, make_list_from_vcf_without_filter
from scripts.HELPERS.paths import create_path_from_GTRD_function
from scripts.HELPERS.paths_for_components import  GTRD_slice_path


interestingSet = {'H1 (embryonic stem cells)',
                  'BGO3 (embryonic stem cells)',
                  'HUES64 (embryonic stem cells)',
                  'CyT49 (embryonic stem cells)',
                  'H9 (embryonic stem cells)',
                  'VAL-3 (embryonic stem cells)',
                  'WA09 (embryonic stem cells)',
                  'embryonic stem cells',
                  'human embryonic stem cells, H1 (WA01)'}

interestingSet = {remove_punctuation(x) for x in interestingSet}
SNP_statistics_dict = {}
with open(GTRD_slice_path, "r") as ml:
    master_list = ml.readlines()
for line in master_list:
    if line[0] == "#":
        continue
    line = line.split("\t")
    if remove_punctuation(line[4]) not in interestingSet:
        continue
    vcf_path = create_path_from_GTRD_function(line, for_what="vcf")
    if not os.path.isfile(vcf_path):
        continue
    list_of_snps = make_list_from_vcf_without_filter(vcf_path)
    for chr, pos, rs_id, ref, alt, ref_counts, alt_counts in list_of_snps:
        try:
            SNP_statistics_dict[(ref_counts, alt_counts)] += 1
        except KeyError:
            SNP_statistics_dict[(ref_counts, alt_counts)] = 1
print(SNP_statistics_dict)
df = pd.DataFrame({'ref': [], 'alt': [], 'count': []})
for ref, alt in SNP_statistics_dict:
    df = df.append(pd.DataFrame({'ref': [ref], 'alt': [alt], 'count': [SNP_statistics_dict[(ref, alt)]]}))
    print(df)
