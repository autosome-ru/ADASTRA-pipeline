import os
import re
import sys
import json

from scripts.HELPERS.paths_for_components import badmaps_path, badmaps_dict_path, correlation_path
from scripts.HELPERS.helpers import Intersection, pack, UnpackBadSegments, get_states


def unpack_snps(line):
    line = line.strip().split('\t')
    return [line[0], int(line[1]), int(line[5]), int(line[6]), line[7]]


def main(file_name):
    if not os.path.isdir(correlation_path):
        try:
            os.mkdir(correlation_path)
        except:
            pass

    with open(badmaps_dict_path, 'r') as file:
        aligns_by_cell_type = json.loads(file.readline().strip())

    modes = []
    for file_name in sorted(os.listdir(badmaps_path)):
        if os.path.isdir(badmaps_path + file_name) and file_name != 'merged_vcfs':
            modes.append(file_name)

    assert os.path.isfile(badmaps_path + 'merged_vcfs/' + file_name)

    name = file_name.split('!')[0]
    lab = file_name.split('!')[1][:-4]

    try:
        aligns = aligns_by_cell_type[file_name[:-4]]  # .tsv
        al_list = [os.path.basename(align) for align in aligns if os.path.isfile(align)]
        datasetsn = len(al_list)
    except KeyError:
        datasetsn = 'nan'
        al_list = []
        print(file_name)

    table_path = badmaps_path + 'merged_vcfs/' + file_name
    for mode in modes:
        if re.match(r'^CAIC@.+@.+$', mode) is not None:
            states = get_states(mode.split('@')[1])
        else:
            states = get_states('')
        if not os.path.isdir(correlation_path + mode + '_tables/'):
            try:
                os.mkdir(correlation_path + mode + '_tables')
            except:
                pass
        badmaps_file_path = badmaps_path + mode + '/' + name + '@' + lab + '.badmap.tsv'
        out_path = correlation_path + mode + '_tables/' + name + '_' + lab.replace('_', '-') + '.tsv'
        print(out_path)

        u = UnpackBadSegments(0)

        with open(table_path, 'r') as table, open(badmaps_file_path, 'r') as ploidy, open(out_path, 'w') as out:
            out.write('#' + str(datasetsn) + '@' + lab + '@' + ','.join(al_list) + '\n')
            for chr, pos, ref, alt, filename, in_intersection, segment_ploidy, segment_id, Qual, segn, sumcov \
                    in Intersection(table, ploidy,
                                    unpack_segments_function=lambda x: u.unpackBADSegments(x, states), unpack_snp_function=unpack_snps,
                                    write_intersect=True, write_segment_args=True):
                if not in_intersection:
                    continue
                out.write(pack([chr, pos, ref, alt, segment_ploidy] + [Qual[x] for x in Qual] + [segn, sumcov] + [filename, segment_id]))


if __name__ == '__main__':
    main(sys.argv[1])
