import sys
import os
from scripts.HELPERS.helpers import ChromPos, pack


def main(peak_path, out_path, p_type):
    with open(out_path, 'w') as o, open(peak_path) as lines:
        for line in lines:
            if line[0] == "#":
                continue
            split_line = line.strip('\n').split('\t')
            split_line[0] = 'chr' + split_line[0]
            if split_line[0] not in ChromPos.chromosomes:
                continue
            if int(split_line[1]) < 0 or int(split_line[2]) < 0:
                continue
            if p_type == "gem":
                split_line[1] = str(max(int(split_line[1]) - 150, 1))
                split_line[2] = str(min(int(split_line[2]) + 150, ChromPos.chromosomes[split_line[0]]))
            o.write(pack(split_line[:3]))


if __name__ == '__main__':
    file = sys.argv[1]
    out = sys.argv[2]
    peak_type = sys.argv[3]
    main(file, out, peak_type)

