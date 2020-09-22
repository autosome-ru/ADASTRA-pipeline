import os
import sys
import json

sys.path.insert(1, "/home/abramov/segmentationValidation/ADASTRA-pipeline")
from scripts.HELPERS.paths_for_components import ploidy_path, ploidy_dict_path, correlation_path
from scripts.HELPERS.helpers import Intersection, pack, unpackBADSegments


def get_states(states_sign):
    if states_sign == '1236':
        states = [1, 2, 3, 6]
    elif states_sign == '12345':
        states = [1, 2, 3, 4, 5]
    elif states_sign == '12345_1.5':
        states = [1, 2, 3, 4, 5, 1.5]
    elif states_sign == '123456':
        states = [1, 2, 3, 4, 5, 6]
    elif states_sign == 'all_but_1.33':
        states = [1, 2, 3, 4, 5, 1.5, 6, 2.5]
    elif states_sign == 'all_but_2.5':
        states = [1, 2, 3, 4, 5, 1.5, 6, 4/3]
    elif states_sign == 'all':
        states = [1, 2, 3, 4, 5, 1.5, 6, 4/3, 2.5]
    else:
        raise ValueError
    return states


def unpack_ploidy_segments(line):
    if line[0] == '#':
        return [''] * 9
    line = line.strip().split('\t')

    return [line[0], int(line[1]), int(line[2]), float(line[3]), int(line[4]),
            int(line[5]), int(line[6]), int(line[7]), int(line[8])]


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

    assert os.path.isfile(ploidy_path + 'merged_vcfs/' + file_name)

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

    table_path = ploidy_path + 'merged_vcfs/' + file_name
    for mode in modes:
        states = get_states(mode.split('@')[1])
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
            for chr, pos, ref, alt, in_intersection, segment_ploidy, Qual, segn, sumcov \
                    in Intersection(table, ploidy,
                                    unpack_segments_function=lambda x: unpackBADSegments(x, [1, 2, 3]), unpack_snp_function=unpack_snps,
                                    write_intersect=True, write_segment_args=True):
                if not in_intersection:
                    continue
                out.write(pack([chr, pos, ref, alt, segment_ploidy] + [Qual[x] for x in Qual] + [segn, sumcov]))
