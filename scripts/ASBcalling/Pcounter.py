import json
import sys
from scipy.stats import binom_test

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import ploidy_dict_path, create_ploidy_path_function
from scripts.HELPERS.helpers import callers_names, unpack, pack, Intersection


def count_p(x, n, p, alternative):
    pv = (binom_test(x, n, p, alternative) + binom_test(x, n, 1 - p, alternative)) / 2
    return pv


def make_reverse_dict(dictionary):
    new_dict = {}
    for key in dictionary:
        paths = dictionary[key]
        for path in paths:
            if path.split("/")[-3] != "CTRL":
                new_dict[path] = key
    return new_dict


if __name__ == '__main__':
    full_path = sys.argv[1]
    
    key = full_path + ".vcf.gz"
    table_annotated = full_path + "_table_annotated.txt"
    output = full_path + "_table_p.txt"
    
    with open(ploidy_dict_path, "r") as read_file:
        d = json.loads(read_file.readline())
        rev_d = make_reverse_dict(d)
    
    ploidy = create_ploidy_path_function(rev_d[key])
    
    print('Now doing {} \n with ploidy file {}'.format(table_annotated, rev_d[key]))
    
    with open(ploidy, 'r') as ploidy_file, open(output, 'w') as out, open(table_annotated, 'r') as table_file:
        
        for chr, pos, ID, ref, alt, ref_c, alt_c, Q, GQ, in_callers, ploidy, dip_qual, lq, rq, seg_c in \
                Intersection(table_file, ploidy_file, write_segment_args=True,
                             unpack_snp_function=lambda x: unpack(x, use_in='Pcounter')):
            if ploidy == '':
                ploidy = 0
                p_ref = 0
                p_alt = 0
            elif float(ploidy) != 0:
                # p_value counting
                p = 1 / (float(ploidy) + 1)
                n = ref_c + alt_c
                p_ref = count_p(ref_c, n, p, 'greater')
                p_alt = count_p(ref_c, n, p, 'less')
            else:
                p_ref = '.'
                p_alt = '.'
            
            out.write(pack([chr, pos, ID, ref, alt, ref_c, alt_c, Q, GQ] +
                           [in_callers[name] for name in callers_names] +
                           [ploidy, dip_qual, lq, rq, seg_c, p_ref, p_alt]))
