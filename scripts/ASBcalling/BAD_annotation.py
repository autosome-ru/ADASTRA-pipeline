import json
import sys

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import create_ploidy_path_function
from scripts.HELPERS.paths_for_components import ploidy_dict_path
from scripts.HELPERS.helpers import callers_names, unpack, pack, Intersection, unpackBADSegments, states


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
    output = full_path + "_table_BADs.txt"

    with open(ploidy_dict_path, "r") as read_file:
        d = json.loads(read_file.readline())
        rev_d = make_reverse_dict(d)

    ploidy_file_name = rev_d[key]

    print('Now doing {} \n with ploidy file {}'.format(table_annotated, ploidy_file_name))

    model = 'CAIC'
    ploidy = create_ploidy_path_function(ploidy_file_name, model)

    with open(ploidy, 'r') as ploidy_file, open(output, 'w') as out, open(table_annotated, 'r') as table_file:
        out.write(pack(['#chr', 'pos', 'ID', 'ref', 'alt', 'ref_read_counts', 'alt_read_counts',
                        'repeat_type'] + callers_names + ['BAD'] + ["Q{:.2f}".format(x) for x in states] +
                       ['SNP_count', 'sum_cover']))

        for chr, pos, ID, ref, alt, ref_c, alt_c, repeat_type, in_callers, \
            in_intersection, BAD, Quals, seg_c, sum_cov in \
                Intersection(table_file, ploidy_file, write_segment_args=True, write_intersect=True,
                             unpack_snp_function=lambda x: unpack(x, use_in='Pcounter'),
                             unpack_segments_function=unpackBADSegments):
            if in_intersection:
                out.write(pack([chr, pos, ID, ref, alt, ref_c, alt_c, repeat_type] +
                               [in_callers[name] for name in callers_names] +
                               [BAD] + [Quals[x] for x in Quals] + [seg_c, sum_cov]))
