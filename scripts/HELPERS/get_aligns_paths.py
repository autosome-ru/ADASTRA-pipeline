import re
import glob
from pathlib import Path
from sys import argv


REGEX_ALIGNS = r"D?ALIGNS\d{4,}"
BASE_PATH_ALIGNS = "/srv/*/egrid/aligns-sorted/*.bam"


def pack(values):
    return '\t'.join(map(str, values)) + '\n'


def find_aligns_paths(aligns):
    paths = {x: None for x in aligns}
    for file in glob.glob(BASE_PATH_ALIGNS):
        if Path(file).is_file():                # If indexed
            align_id = re.search(REGEX_ALIGNS, file)  # and named in db
            if align_id is None:
                continue
            align_id = align_id.group()
            if align_id in aligns:
                paths[align_id] = file
    return paths


def get_files(data):
    aligns_paths = find_aligns_paths([x for x in data['ALIGNS']])
    data['DOWNLOAD_PATH'] = [aligns_paths[x] for x in data['ALIGNS']]
    return data


if __name__ == '__main__':
    with open(argv[1], "r") as infile:
        header = infile.readline().strip('\n').split('\t')
        lines_list = [line.strip('\n').split('\t') for line in infile]
    data_list = dict(zip(header, lines_list))
    data_with_paths = get_files(data_list)
    print('\t'.join(data_with_paths.keys()))
    for i in range(len(data_list['ALIGNS'])):
        print('\t'.join([data_list[x][i] for x in data_list]))


