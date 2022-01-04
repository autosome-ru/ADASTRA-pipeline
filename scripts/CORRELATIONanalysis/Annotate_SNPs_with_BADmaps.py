import os
import re
import sys
import json
from scipy import stats as st
import errno


from scripts.HELPERS.paths_for_components import badmaps_path, badmaps_dict_path
from scripts.HELPERS.helpers import Intersection, pack, UnpackBadSegments, get_states_from_model_name, get_models_list
from scripts.HELPERS.paths import get_badmaps_path_by_validity, get_correlation_path


def unpack_snps(line):
    line = line.strip().split('\t')
    return [line[0], int(line[1]), int(line[5]), int(line[6]), line[7]]


def get_p_value(n, p, x):
    dist = st.binom(n=n, p=p)
    cdf = dist.cdf
    return (cdf(x) - cdf(4) + cdf(n-5) - cdf(n-x-1)) / (cdf(n-5) - cdf(4)) if x < n/2 else 1


def main(file_name, remake=False):
    correlation_path = get_correlation_path()
    with open(badmaps_dict_path, 'r') as file:
        aligns_by_cell_type = json.loads(file.readline().strip())

    modes = get_models_list()

    try:
        assert os.path.isfile(os.path.join(badmaps_path, 'merged_vcfs', file_name))
    except AssertionError:
        print(os.path.join(badmaps_path, 'merged_vcfs', file_name), file_name)
        exit(1)

    name = file_name.split('@')[0]
    lab = os.path.splitext(file_name.split('@')[1])[0]

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
        states = get_states_from_model_name(mode)
        out_dir = os.path.join(correlation_path, mode + '_tables{}'.format('_filtered' if remake else ''))
        if not os.path.isdir(out_dir):
            try:
                os.mkdir(out_dir)
            except OSError as exc:
                if exc.errno != errno.EEXIST:
                    raise
                pass
        badmaps_file_path = os.path.join(get_badmaps_path_by_validity(valid=remake), mode, name + '@' + lab + '.badmap.tsv')
        out_path = os.path.join(out_dir, name + '@' + lab + '.tsv')
        print(out_path)

        u = UnpackBadSegments(0)

        with open(table_path, 'r') as table, open(badmaps_file_path, 'r') as BADmap_file, open(out_path, 'w') as out:
            out.write('#' + str(datasetsn) + '@' + lab + '@' + ','.join(al_list) + '\n')
            for chrom, pos, ref, alt, filename, in_intersection, segment_BAD, segment_snps, segment_snp_ids,\
                    segment_sumcov, segment_id, Qual \
                    in Intersection(table, BADmap_file,
                                    unpack_segments_function=lambda x: u.unpack_bad_segments(x, states),
                                    unpack_snp_function=unpack_snps,
                                    write_intersect=True, write_segment_args=True):
                if not in_intersection:
                    continue
                # p_value = get_p_value(ref + alt, 1 / (segment_BAD + 1), min(ref, alt))
                p_value = 0.5  # FIXME
                out.write(pack([chrom, pos, ref, alt, segment_BAD] +
                               [Qual[x] for x in Qual] + [segment_snp_ids, segment_sumcov] +
                               [filename, segment_id, p_value]))


if __name__ == '__main__':
    main(sys.argv[1])
