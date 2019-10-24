import sys
import gzip
import os.path

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.helpers import pack, make_dict_from_vcf, Intersection, callers_names


# FIXME: change return
def unpack_repeats(line):
    if line[0] == '#':
        return []
    return line.strip().split("\t")


#FIXME
def make_sorted_caller_path(line):
    return line


if __name__ == "__main__":
    exp = dict()
    with gzip.open(sys.argv[1], 'rt') as f:
        make_dict_from_vcf(f, exp)

    sorted_lines = [[chr, pos, ID, REF, ALT, R, A] for ((chr, pos, ID, REF, ALT), (R, A)) in exp.items()]
    sorted_lines = sorted(sorted_lines, key=lambda x: x[1])
    sorted_lines = sorted(sorted_lines, key=lambda x: x[0])

    with open(sys.argv[3], "r") as repeats_file:
        new_arr = []
        for chr, pos, ID, REF, ALT, R, A, in_repeats, repeat_type \
                in Intersection(sorted_lines, repeats_file, write_intersect=True, write_segment_args=True,
                                unpack_segments_function=unpack_repeats):
            if in_repeats and ID == ".":
                continue
            new_arr.append([chr, pos, ID, REF, ALT, R, A, repeat_type])
        sorted_lines = new_arr

    for peak_type in callers_names:
        new_arr = []
        caller_path = make_sorted_caller_path(callers_names)
        if os.path.isfile(caller_path):
            peak_file = open(caller_path, "r")
        else:
            peak_file = []
        for chr, pos, ID, REF, ALT, R, A, repeat_type, *in_peaks in Intersection(sorted_lines, peak_file,
                                                                                 write_intersect=True):
            new_arr.append([chr, pos, ID, REF, ALT, R, A, repeat_type] + in_peaks)
        sorted_lines = new_arr

    with open(sys.argv[2], "w") as out:
        for split_line in sorted_lines:
            out.write(pack(split_line))
