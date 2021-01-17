import re
import glob
from pathlib import Path
from sys import argv

if argv[1] == 'chip':
  REGEX_ALIGNS = r"ALIGNS\d{4,}"
elif argv[1] == 'dnase':
  REGEX_ALIGNS = r"DALIGNS\d{4,}"
elif argv[1] == 'atac':
    REGEX_ALIGNS = r"AALIGNS\d{4,}"
elif argv[1] == 'faire':
    REGEX_ALIGNS = r"FALIGNS\d{4,}"
else:
    raise ValueError('Must be chip or dnase')
BASE_PATH_ALIGNS = "/srv/*/egrid/aligns-sorted/*.bam"


def pack(values):
    return '\t'.join(map(str, values)) + '\n'


def find_aligns_paths(aligns):
    paths = {x: [] for x in aligns}
    for file in glob.glob(BASE_PATH_ALIGNS):
        if Path(file).is_file():                # If indexed
            align_id = re.search(REGEX_ALIGNS, file)  # and named in db
            if align_id is None:
                continue
            align_id = align_id.group()
            if align_id in aligns:
                paths[align_id].append(file)
    for key, value in paths.values():
        paths[key] = ' '.join(value) if value else None
    return paths


def get_files(data):
    aligns_paths = find_aligns_paths(data['ALIGNS'])
    data['DOWNLOAD_PATH'] = [aligns_paths[x] for x in data['ALIGNS']]
    return data


if __name__ == '__main__':
    data_list = {}
    with open(argv[2], "r") as infile:
        header = infile.readline().strip('\n').split('\t')
        lines_list = [line.strip('\n').split('\t') for line in infile if line.strip()]

    for index, item in enumerate(header):
        data_list[item] = [x[index] for x in lines_list]
    data_list = get_files(data_list)
    header.append('DOWNLOAD_PATH')
    print('\t'.join(header))
    for i in range(len(data_list['ALIGNS'])):
        print('\t'.join(map(str, [data_list[x][i] for x in header])))