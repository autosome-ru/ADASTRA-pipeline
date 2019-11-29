import os
import sys
import json

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import ploidy_dict_path, ploidy_path, correlation_path
from scripts.HELPERS.helpers import Intersection, pack, read_synonims, corrected_BADs_dict


def unpack_ploidy_segments(line):
    if line[0] == '#':
        return [''] * 6
    line = line.strip().split('\t')

    return [line[0], int(line[1]), int(line[2]), float(line[3]), int(line[4]), int(line[7]), int(line[8])]


def unpack_snps(line):
    line = line.strip().split('\t')
    return [line[0], int(line[1]), int(line[5]), int(line[6])]


if __name__ == '__main__':
    if not os.path.isdir(correlation_path):
        try:
            os.mkdir(correlation_path)
        except:
            pass

    with open(ploidy_dict_path, 'r') as file:
        aligns_by_cell_type = json.loads(file.readline().strip())

    modes = []
    for file_name in sorted(os.listdir(ploidy_path)):
        if os.path.isdir(ploidy_path + file_name) and file_name != 'merged_vcfs':
            modes.append(file_name)

    file_name = sys.argv[1]

    assert os.path.isfile(ploidy_path + file_name)

    name = file_name.split('!')[0]
    lab = file_name.split('!')[1][:-4]

    try:
        aligns = aligns_by_cell_type[file_name[:-4]]  # .tsv
        al_list = [align[29:-7] for align in aligns if os.path.isfile(align)]
        datasetsn = len(al_list)
    except KeyError:
        datasetsn = 'nan'
        al_list = []
        print(file_name)

    names, _ = read_synonims()
    if name in names:
        table_path = ploidy_path + 'merged_vcfs/' + file_name
        for mode in modes:
            if not os.path.isdir(correlation_path + mode + '_tables/'):
                try:
                    os.mkdir(correlation_path + mode + '_tables')
                except:
                    pass
            ploidy_file_path = ploidy_path + mode + '/' + name + '!' + lab + '_ploidy.tsv'
            out_path = correlation_path + mode + '_tables/' + name + '_' + lab.replace('_', '-') + '.tsv'
            print(out_path)

            with open(table_path, 'r') as table, open(ploidy_file_path, 'r') as ploidy, open(out_path, 'w') as out:
                out.write('#' + str(datasetsn) + '!' + lab + '!' + '>'.join(al_list) + '\n')
                for chr, pos, ref, alt, in_intersection, segment_ploidy, qual, segn, sumcov \
                        in Intersection(table, ploidy,
                                        unpack_segments_function=unpack_ploidy_segments, unpack_snp_function=unpack_snps,
                                        write_intersect=True, write_segment_args=True):
                    if not in_intersection:
                        continue
                    out.write(pack([chr, pos, ref, alt, corrected_BADs_dict.get(segment_ploidy, segment_ploidy),
                                    qual, segn, sumcov]))
