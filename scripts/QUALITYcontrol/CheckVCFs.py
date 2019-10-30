import gzip
import os.path
import sys

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import parameters_path
from scripts.HELPERS.helpers import ChromPos, make_list_for_VCFs

out_path = parameters_path + "BadVCFs.tsv"
chrs = dict(zip(['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY'], [i for i in range(24)]))


def check_vcf(path, missing_chromosomes_threshold=2):
    a = [False] * 24
    is_bad_vcf = False
    if not os.path.isfile(path):
        return False
    with gzip.open(path, 'rt') as vcf:
        for line in vcf:
            if line[0] == '#':
                continue
            line = line.strip().split("\t")
            chr = line[0]
            if chr not in ChromPos.chrs:
                is_bad_vcf = True
                break
            if not a[chrs[chr]]:
                a[chrs[chr]] = True

    number_of_bad_chromosomes = 0
    print(path)
    if is_bad_vcf:
        return is_bad_vcf
    for chr_index in reversed(a):
        print(chr_index)
        if not chr_index:
            number_of_bad_chromosomes += 1
        else:
            if number_of_bad_chromosomes > missing_chromosomes_threshold:
                is_bad_vcf = True
            break
    return is_bad_vcf


def check_chromosome(path, chromosome_name):
    if not os.path.isfile(path):
        return False
    with gzip.open(path, 'rt') as vcf:
        for line in vcf:
            if line[0] == '#':
                continue
            line = line.strip().split("\t")
            chr = line[0]
            if chr not in ChromPos.chrs:
                print(chr)
            if chr == chromosome_name:
                print(line)
                return True
    return False


if __name__ == "__main__":
    d = make_list_for_VCFs(condition_function=lambda x: check_chromosome(x, chromosome_name="chr11"))
    print("{}/{} vcfs have 11 chromosome".format(sum(x[1] for x in d.items()), len(d)))
    d = make_list_for_VCFs(condition_function=lambda x: check_chromosome(x, chromosome_name="chr1"))
    print("{}/{} vcfs have 1 chromosome".format(sum(x[1] for x in d.items()), len(d)))
