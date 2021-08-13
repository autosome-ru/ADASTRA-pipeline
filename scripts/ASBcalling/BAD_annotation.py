import json
import os
import sys
from scripts.HELPERS.paths import create_badmaps_path_function, get_ending, get_merged_badmaps_dict_path
from scripts.HELPERS.helpers import callers_names, unpack, pack, Intersection, UnpackBadSegments, segmentation_states


def make_reverse_dict(dictionary):
    new_dict = {}
    for key in dictionary:
        paths = dictionary[key]
        for path in paths:
            new_dict[path] = key
    return new_dict


def main(key, remade=True):
    table_annotated = key + get_ending("annotation")
    output = key + get_ending("BAD")

    with open(get_merged_badmaps_dict_path(remade=remade), "r") as read_file:
        d = json.load(read_file)
        rev_d = make_reverse_dict(d)

    badmap_file_name = rev_d[key]

    print('Now doing {} \n with BAD map file {}'.format(table_annotated, badmap_file_name))
    badmap_file_path = create_badmaps_path_function(badmap_file_name, valid=True)
    with open(badmap_file_path, 'r') as badmap_file, open(output, 'w') as out, open(table_annotated, 'r') as table_file:
        out.write(pack(['#chr', 'pos', 'ID', 'ref', 'alt', 'ref_read_counts', 'alt_read_counts',
                        'repeat_type'] + callers_names + ['BAD'] + ["Q{:.2f}".format(x) for x in segmentation_states] +
                       ['SNP_count', 'sum_cover']))

        u = UnpackBadSegments(None)
        for chr, pos, ID, ref, alt, ref_c, alt_c, repeat_type, in_callers, \
            in_intersection, segment_BAD, segment_snps, segment_snp_ids,\
                    segment_sumcov, Qual in \
                Intersection(table_file, badmap_file, write_segment_args=True, write_intersect=True,
                             unpack_snp_function=lambda x: unpack(x, use_in='Pcounter'),
                             unpack_segments_function=lambda x: u.unpack_bad_segments(x, segmentation_states)):
            if in_intersection and ID.startswith('rs') and segment_BAD != 6.0:
                out.write(pack([chr, pos, ID, ref, alt, ref_c, alt_c, repeat_type] +
                               [in_callers[name] for name in callers_names] +
                               [segment_BAD] + [Qual[x] for x in Qual] + [segment_snp_ids, segment_sumcov]))


if __name__ == '__main__':
    main(sys.argv[1])
