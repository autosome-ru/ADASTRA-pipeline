import sys
import gzip

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.helpers import make_dict_from_vcf


# TODO FIX Annotate.py

if __name__ == "__main__":
    vcf = gzip.open(sys.argv[1], 'rt')

    exp = dict()
    make_dict_from_vcf(vcf, exp)

    exp_keys = list(exp.keys())
    exp_keys = sorted(exp_keys, key=lambda x: x[1])
    exp_keys = sorted(exp_keys, key=lambda x: x[0])

    # TODO: check if in repeat regions
    out = sys.argv[2]
