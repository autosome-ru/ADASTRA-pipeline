import sys
import gzip
import os.path
from scripts.HELPERS.helpers import pack, make_dict_from_vcf, Intersection, callers_names
from scripts.HELPERS.paths_for_components import repeats_path

from scripts.HELPERS.paths import get_ending


def make_sorted_caller_path(base_path, peak_type):
    return '{}_{}.bed.sorted'.format(base_path, peak_type)


def main(base_path):
    exp = dict()
    with gzip.open(base_path + get_ending('vcf'), 'rt') as f:
        make_dict_from_vcf(f, exp)
    sorted_lines = [[chr, pos, ID, REF, ALT, R, A] for ((chr, pos, ID, REF, ALT), (R, A)) in exp.items()]
    sorted_lines = sorted(sorted_lines, key=lambda x: x[1])
    sorted_lines = sorted(sorted_lines, key=lambda x: x[0])
    with open(repeats_path, "r") as repeats_buffer:
        new_arr = []
        for chr, pos, ID, REF, ALT, R, A, in_repeats, repeat_type \
                in Intersection(sorted_lines, repeats_buffer, write_intersect=True, write_segment_args=True):
            if in_repeats and ID == ".":
                continue
            new_arr.append([chr, pos, ID, REF, ALT, R, A, repeat_type])
        sorted_lines = new_arr

    for peak_type in callers_names:
        new_arr = []
        caller_path = make_sorted_caller_path(base_path, peak_type)
        if os.path.isfile(caller_path):
            peak_file = open(caller_path, "r")
        else:
            peak_file = []
        for chr, pos, ID, REF, ALT, R, A, repeat_type, *in_peaks in Intersection(sorted_lines, peak_file,
                                                                                 write_intersect=True):
            new_arr.append([chr, pos, ID, REF, ALT, R, A, repeat_type] + in_peaks)
        sorted_lines = new_arr
    table_annotated_path = base_path + get_ending('annotation')
    with open(table_annotated_path, "w") as out:
        out.write(pack(['#chr', 'pos', 'ID', 'ref', 'alt', 'ref_read_counts', 'alt_read_counts', 'repeat_type'] +
                       callers_names))
        for split_line in sorted_lines:
            out.write(pack(split_line))


if __name__ == '__main__':
    main(sys.argv[1])
