import re
import glob
from pathlib import Path
from sys import argv

REGEX_PEAKS = r"PEAKS\d{4,}"
REGEX_ALIGNS = r"D?ALIGNS\d{4,}"

BASE_PATH_ALIGNS = "/srv/*/egrid/aligns-sorted/*.bam"


def extract_aligns(aligns):
    paths = {x: "None" for x in aligns}
    for file in glob.glob(BASE_PATH_ALIGNS):
        if Path(file).is_file():                # If indexed
            id = re.search(REGEX_ALIGNS, file)  # and named in db
            if id is None:
                continue
            id = id.group()
            if id in paths:
                paths[id] = file
    return paths


def get_files(lines):
    aligns = {}
    for line in lines:
        if line[0] == "#":
            continue
        line = line.strip().split("\t")
        aligns[line[6]] = line
    aligns_paths = extract_aligns(aligns)

    for align_name in aligns:
        print("\t".join(aligns[align_name] + [aligns_paths[align_name]]))


with open(argv[1], "r") as infile:
    get_files(infile.readlines())
