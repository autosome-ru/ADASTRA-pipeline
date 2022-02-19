import gzip
import os
import sys
import json

from scripts.HELPERS.paths import get_ending, create_merged_vcf_path_function, get_new_badmaps_dict_path
from scripts.HELPERS.paths_for_components import badmaps_dict_path
from scripts.HELPERS.helpers import make_dict_from_vcf, make_list_from_vcf, default_model


def merge_vcfs_add_counts(out_file_name, in_files):
    vcf_dict = dict()
    for file in in_files:
        with gzip.open(file, 'rt') as vcf:
            make_dict_from_vcf(vcf, vcf_dict)

    vcf_keys = list(vcf_dict.keys())
    vcf_keys.sort(key=lambda cords: cords[1])
    vcf_keys.sort(key=lambda cords: cords[0])

    with open(out_file_name, 'w') as out:
        for (chr, pos, ID, REF, ALT) in vcf_keys:
            (R, A) = vcf_dict[(chr, pos, ID, REF, ALT)]
            out.write('\t'.join(map(str, [chr, pos, ID, REF, ALT, R, A])) + '\n')


def merge_vcfs_as_independent_snps(out_file_name, in_files):
    vcf_list = []
    for file in in_files:
        with gzip.open(file, 'rt') as vcf:
            vcf_list += make_list_from_vcf(vcf, file_name=os.path.basename(file))

    vcf_list.sort(key=lambda cords: cords[1])
    vcf_list.sort(key=lambda cords: cords[0])

    with open(out_file_name, 'w') as out:
        for (chr, pos, ID, REF, ALT, R, A, filename) in vcf_list:
            out.write('\t'.join(map(str, [chr, pos, ID, REF, ALT, R, A, filename])) + '\n')


def main(key, remake=False):
    with open(get_new_badmaps_dict_path(default_model) if remake else badmaps_dict_path, 'r') as read_file:
        d = json.load(read_file)
    mode = 'independent'

    paths_list = []
    for path in d[key]:
        if os.path.isfile(path + get_ending("vcf")):
            paths_list.append(path + get_ending("vcf"))
    out_file = create_merged_vcf_path_function(key)

    if mode == 'independent':
        merge_vcfs_as_independent_snps(out_file, paths_list)
    elif mode == 'add':
        merge_vcfs_add_counts(out_file, paths_list)
    else:
        raise ValueError(mode)


if __name__ == '__main__':
    main(sys.argv[1])
