import io
import sys
import zipfile

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
            o.write(line)
            continue
        split_line = line.split()
        if int(split_line[1]) < 0 or int(split_line[2]) < 0:
            continue
        o.write(line)
