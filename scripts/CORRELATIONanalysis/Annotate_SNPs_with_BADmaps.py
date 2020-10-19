import os
import re
import sys
import json
from scipy import stats as st
import numpy as np
import pandas as pd

from scripts.HELPERS.paths_for_components import badmaps_path, badmaps_dict_path, correlation_path
from scripts.HELPERS.helpers import Intersection, pack, UnpackBadSegments, get_states


def unpack_snps(line):
    line = line.strip().split('\t')
    return [line[0], int(line[1]), int(line[5]), int(line[6]), line[7]]


def get_p_value(n, p, x):
    dist = st.binom(n=n, p=p)
    cdf = dist.cdf
    return (cdf(x) - cdf(4) + cdf(n-5) - cdf(n-x-1)) / (cdf(n-5) - cdf(4)) if x < n/2 else 1


def main(file_name):
    if not os.path.isdir(correlation_path):
        try:
            os.mkdir(correlation_path)
        except:
            pass

    with open(badmaps_dict_path, 'r') as file:
        aligns_by_cell_type = json.loads(file.readline().strip())

    modes = []
    for dir_name in sorted(os.listdir(badmaps_path)):
        if os.path.isdir(os.path.join(badmaps_path, dir_name)) and dir_name not in ('merged_vcfs', 'valid_BADmaps'):
            modes.append(dir_name)

    try:
        assert os.path.isfile(os.path.join(badmaps_path, 'merged_vcfs', file_name))
    except AssertionError:
        print(os.path.join(badmaps_path, 'merged_vcfs', file_name), file_name)
        exit(1)

    name = file_name.split('@')[0]
    lab = file_name.split('@')[1][:-4]  # .tsv

    try:
        aligns = aligns_by_cell_type[file_name[:-4]]  # .tsv
        al_list = [os.path.basename(align) for align in aligns if os.path.isfile(align)]
        datasetsn = len(al_list)
    except KeyError:
        datasetsn = 'nan'
        al_list = []
        print(file_name)

    table_path = os.path.join(badmaps_path, 'merged_vcfs', file_name)
    for mode in modes:
        if re.match(r'^CAIC@.+@.+$', mode) is not None:
            states = get_states(mode.split('@')[1])
        else:
            states = get_states('')
        if not os.path.isdir(os.path.join(correlation_path, mode + '_tables')):
            try:
                os.mkdir(os.path.join(correlation_path, mode + '_tables'))
            except:
                pass
        badmaps_file_path = os.path.join(badmaps_path, mode, name + '@' + lab + '.badmap.tsv')
        out_path = os.path.join(correlation_path, mode + '_tables', name + '@' + lab + '.tsv')
        print(out_path)

        u = UnpackBadSegments(0)

        with open(table_path, 'r') as table, open(badmaps_file_path, 'r') as BADmap_file, open(out_path, 'w') as out:
            out.write('#' + str(datasetsn) + '@' + lab + '@' + ','.join(al_list) + '\n')
            for chr, pos, ref, alt, filename, in_intersection, segment_BAD, segment_id, Qual, segn, sumcov \
                    in Intersection(table, BADmap_file,
                                    unpack_segments_function=lambda x: u.unpackBADSegments(x, states), unpack_snp_function=unpack_snps,
                                    write_intersect=True, write_segment_args=True):
                if not in_intersection:
                    continue
                p_value = get_p_value(ref + alt, 1 / (segment_BAD + 1), min(ref, alt))
                out.write(pack([chr, pos, ref, alt, segment_BAD] + [Qual[x] for x in Qual] + [segn, sumcov] + [filename, segment_id, p_value]))


if __name__ == '__main__':
    main(sys.argv[1])
