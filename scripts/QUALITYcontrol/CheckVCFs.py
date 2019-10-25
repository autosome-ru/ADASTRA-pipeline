import gzip
import os.path
import sys

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import create_path_from_GTRD_function, make_black_list, GTRD_slice_path, parameters_path, \
    create_line_for_snp_calling


out_path = parameters_path + "BadVCFs.tsv"

chrs = dict(zip(['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY'], [i for i in range(24)]))


def check_vcf(path, missing_chromosomes_threshold=3):
    a = [False] * 24
    with gzip.open(path, 'rt') as vcf:
        for line in vcf:
            if line[0] == '#':
                continue
            line = line.split()
            chr = line[0]
            if not a[chrs[chr]]:
                a[chrs[chr]] = True

    is_bad_vcf = False
    number_of_bad_chromosomes = 0
    print(path)
    for chr_index in reversed(a):
        print(chr_index)
        if not chr_index:
            number_of_bad_chromosomes += 1
        else:
            if number_of_bad_chromosomes > missing_chromosomes_threshold:
                is_bad_vcf = True
            break
    return is_bad_vcf


black_list = make_black_list()
counted_controls = set()
if __name__ == "__main__":
    with open(GTRD_slice_path, "r") as master_list, open(out_path, "w") as out:
        for line in master_list:
            if line[0] == "#":
                continue
            split_line = line.strip().split("\t")
            if split_line[0] not in black_list:
                vcf_path = create_path_from_GTRD_function(split_line, for_what="vcf")
                if os.path.isfile(vcf_path):
                    if check_vcf(vcf_path):
                        out.write(create_line_for_snp_calling(split_line))
            if len(split_line) > 10 and split_line[10] not in black_list:
                vcf_path = create_path_from_GTRD_function(split_line, for_what="vcf", ctrl=True)
                if vcf_path in counted_controls:
                    continue
                counted_controls.add(vcf_path)
                if os.path.isfile(vcf_path):
                    if check_vcf(vcf_path):
                        out.write(create_line_for_snp_calling(split_line, is_ctrl=True))
