import json
import sys
from scripts.HELPERS.paths import create_badmaps_path_function, get_ending
from scripts.HELPERS.paths_for_components import badmaps_dict_path
from scripts.HELPERS.helpers import callers_names, unpack, pack, Intersection, UnpackBadSegments, states


def make_reverse_dict(dictionary):
    new_dict = {}
    for key in dictionary:
        paths = dictionary[key]
        for path in paths:
            new_dict[path] = key
    return new_dict


def main(key):
    table_annotated = key + get_ending("annotation")
    output = key + get_ending("BAD")

    with open(badmaps_dict_path, "r") as read_file:
        d = json.loads(read_file.readline())
        rev_d = make_reverse_dict(d)

    badmap_file_name = rev_d[key]

    print('Now doing {} \n with BAD map file {}'.format(table_annotated, badmap_file_name))
    badmap_file_path = create_badmaps_path_function(badmap_file_name, valid=True)

    with open(badmap_file_path, 'r') as badmap_file, open(output, 'w') as out, open(table_annotated, 'r') as table_file:
        out.write(pack(['#chr', 'pos', 'ID', 'ref', 'alt', 'ref_read_counts', 'alt_read_counts',
                        'repeat_type'] + callers_names + ['BAD'] + ["Q{:.2f}".format(x) for x in states] +
                       ['SNP_count', 'sum_cover']))

        u = UnpackBadSegments(None)
        for chr, pos, ID, ref, alt, ref_c, alt_c, repeat_type, in_callers, \
            in_intersection, BAD, Quals, seg_c, sum_cov in \
                Intersection(table_file, badmap_file, write_segment_args=True, write_intersect=True,
                             unpack_snp_function=lambda x: unpack(x, use_in='Pcounter'),
                             unpack_segments_function=lambda x: u.unpackBADSegments(x, states)):
            if in_intersection and ID.startswith('rs'):
                out.write(pack([chr, pos, ID, ref, alt, ref_c, alt_c, repeat_type] +
                               [in_callers[name] for name in callers_names] +
                               [BAD] + [Quals[x] for x in Quals] + [seg_c, sum_cov]))


if __name__ == '__main__':
    main(sys.argv[1])
