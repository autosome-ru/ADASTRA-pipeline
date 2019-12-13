import json
import sys
from scipy.stats import binom_test, binom
import pandas as pd


sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import ploidy_dict_path, correlation_path, create_ploidy_path_function
from scripts.HELPERS.helpers import callers_names, unpack, pack, Intersection, sBAD_dict


def count_p(x, n, p, alternative, cut_off=2):
    if x <= cut_off or x >= n - cut_off:
        raise ValueError('Read-counts {} must be greater than cut off value {}'.format(x, cut_off))
    tail = binom_test(cut_off, n, p, alternative='less') + binom_test(n - cut_off, n, p, alternative='greater')

    pv = (binom_test(x, n, p, alternative) + binom_test(x, n, 1 - p, alternative)) / 2
    pv_corrected = 0.5 - 1 / (1 - tail) * 0.5 + pv / (1 - tail)
    pv_balanced = pv_corrected - (binom(n, p).pmf(x) + binom(n, 1-p).pmf(x)) / 4

    return pv, pv_corrected, pv_balanced


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

    ploidy_file_name = rev_d[key]

    print('Now doing {} \n with ploidy file {}'.format(table_annotated, ploidy_file_name))

    cor_stats = pd.read_table(correlation_path + 'cor_stats_test.tsv')
    names = cor_stats['#cell_line'] + '!' + cor_stats['cells']
    cor_stats['names'] = names

    sum_cov = int(cor_stats[cor_stats['names'] == ploidy_file_name]['sum_cov'])

    # if sum_cov < 15000000:
    #     model = 'CAIC'
    # else:
    #     model = 'SQRT'
    model = 'CAIC'
    ploidy = create_ploidy_path_function(ploidy_file_name, model)

    with open(ploidy, 'r') as ploidy_file, open(output, 'w') as out, open(table_annotated, 'r') as table_file:
        out.write(pack(['#chr', 'pos', 'ID', 'ref', 'alt', 'ref_read_counts', 'alt_read_counts',
                        'repeat_type'] + callers_names + ['BAD', 'Q1', 'left_qual', 'right_qual', 'SNP_count',
                                                          'sum_cover',
                                                          'p_value_ref', 'p_value_alt',
                                                          'p_value_ref_corrected', 'p_value_alt_corrected',
                                                          'p_value_ref_balanced', 'p_value_alt_balanced',
                                                          ]))

        for chr, pos, ID, ref, alt, ref_c, alt_c, repeat_type, in_callers, \
            in_intersection, BAD, dip_qual, lq, rq, seg_c, sum_cov in \
                Intersection(table_file, ploidy_file, write_segment_args=True, write_intersect=True,
                             unpack_snp_function=lambda x: unpack(x, use_in='Pcounter')):
            if in_intersection:
                #  p_value counting
                p = 1 / (float(BAD) + 1)
                n = ref_c + alt_c
                p_ref, p_ref_cor, p_ref_bal = count_p(ref_c, n, p, 'greater')
                p_alt, p_alt_cor, p_alt_bal = count_p(ref_c, n, p, 'less')
            else:
                BAD = 0
                p_ref, p_ref_cor, p_ref_bal = '.', '.', '.'
                p_alt, p_alt_cor, p_alt_bal = '.', '.', '.'

            out.write(pack([chr, pos, ID, ref, alt, ref_c, alt_c, repeat_type] +
                           [in_callers[name] for name in callers_names] +
                           [BAD, dip_qual, lq, rq, seg_c, sum_cov,
                            p_ref, p_alt,
                            p_ref_cor, p_alt_cor,
                            p_ref_bal, p_alt_bal,
                            ]))
