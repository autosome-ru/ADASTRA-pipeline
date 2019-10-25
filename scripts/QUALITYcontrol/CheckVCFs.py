import sys
import gzip

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.helpers import make_dict_from_vcf

if __name__ == "__main__":
    a = [False for i in range(24)]
    chrs = dict(zip(['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY'], [i for i in range(24)]))
    with gzip.open(sys.argv[1], 'rt') as vcf:
        for line in vcf:
            if line[0] == '#':
                continue
            line = line.split()
            chr = line[0]
            if not a[chrs[chr]]:
                a[chrs[chr]] = True
#TODO: make it works

