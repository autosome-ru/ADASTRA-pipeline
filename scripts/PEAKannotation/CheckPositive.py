import io
import sys
import zipfile

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.helpers import ChromPos

file = sys.argv[1]
out = sys.argv[2]

name = ".".join(file.split(".")[-3:-1]).split("/")[-1]
with zipfile.ZipFile(file, "r") as archive:
    f = archive.open(name, "r")
    f = io.TextIOWrapper(f)
    lines = f.readlines()

with open(out, 'w') as o:
    for line in lines:
        if line[0] == "#":
            continue
        split_line = line.split()
        split_line[0] = 'chr' + split_line[0]
        if split_line[0] not in ChromPos.chrs:
            continue
        if int(split_line[1]) < 0 or int(split_line[2]) < 0:
            continue
        o.write('\t'.join(split_line[:3]) + '\n')
