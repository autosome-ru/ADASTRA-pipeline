import json
import os
import re
import sys
from scipy.stats import kendalltau
import numpy as np
import pandas as pd
import errno
from scripts.HELPERS.helpers import CorrelationReader, Intersection, pack, read_synonims, ChromPos, get_states, \
    test_percentiles_list, cover_procentiles_list
from scripts.HELPERS.paths_for_components import cgh_path, cosmic_path, badmaps_path
from scripts.HELPERS.paths import get_correlation_path, get_heatmap_data_path, get_badmaps_path_by_validity

correlation_path = get_correlation_path()
heatmap_data_path = get_heatmap_data_path()
valid_badmaps_path = get_badmaps_path_by_validity(True)


def get_name_by_dir(dir_name, naive_modes):
    if dir_name in naive_modes:
        return dir_name
    return dir_name[:dir_name.rfind('_')].split('/')[-1]


def count_cosmic_segments(cosmic_names, cell_line_name):
    with open(cosmic_path, 'r') as cosmic_file:
        return sum(True for line in cosmic_file if unpack_cosmic_segments(line,
                                                                          cell_line_name=cell_line_name,
                                                                          cosmic_names=cosmic_names))


def unpack_cosmic_segments(line, mode='normal', cosmic_names=None, cell_line_name=''):
    if cosmic_names is None:
        cosmic_names = {}
    if line[0] == '#':
        return [''] * len(line.strip().split('\t'))
    line = line.strip().split('\t')
    if cell_line_name not in cosmic_names or line[0] != cosmic_names.get(cell_line_name):
        return []
    if int(line[4]) == 0:
        return []

    if mode == 'normal':
        value = int(line[5]) / int(line[4]) - 1
    elif mode == 'total':
        value = int(line[5])
    else:
        raise ValueError(mode)

    return [line[1], int(line[2]), int(line[3]), value]


def correlation_with_cosmic(SNP_objects, mode, method='normal',
                            heatmap_data_file=None, cell_line_name='',
                            cosmic_names=None):
    if cosmic_names is None:
        cosmic_names = {}
    heatmap = None if heatmap_data_file is None else open(heatmap_data_file, 'w')
    cosmic_segments = []
    with open(cosmic_path, 'r') as cosmic_file:
        for line in cosmic_file:
            cosmic_cell_line_segments = unpack_cosmic_segments(line, mode=mode,
                                                               cell_line_name=cell_line_name,
                                                               cosmic_names=cosmic_names)
            if cosmic_cell_line_segments:
                cosmic_segments.append(cosmic_cell_line_segments)
    snp_BAD_list = []
    cosmic_BAD_list = []
    if method == 'normal':
        for chromosome, pos, snp_BAD, quals, in_intersect, cosmic_BAD \
                in Intersection(SNP_objects, cosmic_segments, write_intersect=True,
                                write_segment_args=True):
            if not in_intersect:
                continue
            snp_BAD_list.append(snp_BAD)
            cosmic_BAD_list.append(cosmic_BAD)

            if heatmap is not None:
                heatmap.write(pack([chromosome, pos, snp_BAD, cosmic_BAD]))
        if heatmap is not None:
            heatmap.close()

        if len(snp_BAD_list) != 0:
            kt = kendalltau(snp_BAD_list, cosmic_BAD_list)[0]
            if kt == 'nan':
                return 'NaN'
            return kt
        return 'NaN'
    elif method == 'cover':
        for chromosome, pos, cov, snp_BAD, quals, in_intersect, cosmic_BAD \
                in Intersection(SNP_objects, cosmic_segments, write_intersect=True,
                                write_segment_args=True):
            if not in_intersect:
                continue
            snp_BAD_list.append(snp_BAD)
            cosmic_BAD_list.append(cosmic_BAD)

            if heatmap is not None:
                heatmap.write(pack([chromosome, pos, cov, snp_BAD, cosmic_BAD] +
                                   [quals[x] for x in quals]))
        if heatmap is not None:
            heatmap.close()

        if len(snp_BAD_list) != 0:
            kt = kendalltau(snp_BAD_list, cosmic_BAD_list)[0]
            if kt == 'nan':
                return 'NaN'
            return kt
        return 'NaN'


def find_nearest_probe_to_SNP(SNP_objects, CGH_objects):
    nearest_probes = []
    if not CGH_objects:
        return []
    i = 0
    for SNP in SNP_objects:
        SNP = [ChromPos(SNP[0], SNP[1])] + SNP[2:]
        current_distance = SNP[0].distance(ChromPos(CGH_objects[i][0], CGH_objects[i][1]))
        while SNP[0].distance(ChromPos(CGH_objects[i + 1][0], CGH_objects[i + 1][1])) <= current_distance:
            current_distance = SNP[0].distance(ChromPos(CGH_objects[i + 1][0], CGH_objects[i + 1][1]))
            i += 1
        nearest_probes.append(CGH_objects[i])
    return nearest_probes


def get_p_value_quantiles(percentiles_list, df):
    return list(np.quantile(df['p_value'], [x / 100 for x in percentiles_list]))


def get_p_value_tails_weights(percentiles_list, df):
    total = len(df.index)
    return list(len(df[df['p_value'] <= p / 100].index) / total for p in percentiles_list)


def filter_segments_or_datasets(snps_path, states, percentiles_list):
    with open(snps_path, 'r') as out:
        if not out.readline():
            return ['NaN'] * (len(percentiles_list) + len(cover_procentiles_list)), '{}'
    out_table = pd.read_table(snps_path, header=None, comment='#')
    out_table.columns = ['chr', 'pos', 'ref', 'alt', 'BAD'] + ['Q{:.2f}'.format(BAD) for BAD in states] + ['snps_n',
                                                                                                           'sumcov',
                                                                                                           'dataset',
                                                                                                           'seg_id',
                                                                                                           'p_value']
    quals = get_p_value_quantiles(test_percentiles_list, out_table)
    p_tails = get_p_value_tails_weights(test_percentiles_list, out_table)
    cover_list = out_table['ref'] + out_table['alt']
    q_var = [np.quantile(cover_list, q/100) for q in cover_procentiles_list]

    datasets_info = {}
    for dataset in out_table['dataset'].unique():
        table = out_table[out_table['dataset'] == dataset]
        quals = get_p_value_quantiles(test_percentiles_list, table)
        cover_list = table['ref'] + table['alt']
        total = len(table.index)
        datasets_info[dataset] = {
            'snps': len(table.index),
            'cover': sum(cover_list),
            'quals': {'Q{}'.format(pr): q for pr, q in zip(test_percentiles_list, quals)},
            'p_tails': {'P{}'.format(p): len(table[table['p_value'] <= p / 100].index) / total for p in percentiles_list},
            'cover_quals': {'CQ{}'.format(q): np.quantile(cover_list, q/100) for q in cover_procentiles_list},
        }

    return quals + p_tails + q_var, datasets_info


def main(file_name, remake=False):
    print(file_name)

    out_path = os.path.join(correlation_path, file_name + '.thread')

    snp_dirs = []
    naive_modes = ['naive']

    for f_name in sorted(os.listdir(correlation_path)):
        snp_dir = os.path.join(correlation_path, f_name)
        if f_name.endswith('_tables') and os.path.isdir(snp_dir):
            snp_dirs.append(snp_dir)

    reader = CorrelationReader()
    reader.CGH_path = cgh_path

    cosmic_names, cgh_names = read_synonims()

    with open(out_path, 'w') as out:

        # if file_name != 'HCT-116_colon_carcinoma_19.tsv': continue

        corr_to_objects = {}
        segment_numbers = {}
        quality_scores = {}
        datasets_info = {}
        # print('reading COSMIC')
        cell_line_name = file_name[:file_name.rfind('@')]

        for snp_dir in snp_dirs:
            model = get_name_by_dir(snp_dir, naive_modes)

            if re.match(r'^CAIC@.+@.+$', model) is not None:
                print(model.split('@')[1])
                states = get_states(model.split('@')[1])
            else:
                states = get_states('')
            reader.states = states

            reader.SNP_path = os.path.join(snp_dir, file_name)

            if remake:
                new_dir = snp_dir[:-1] + '_filtered' if snp_dir.endswith('/') else snp_dir + '_filtered'
                if not os.path.isdir(new_dir):
                    os.mkdir(new_dir)
                new_path = os.path.join(new_dir, file_name)
                reader.SNP_path = new_path

            quality_scores[model], datasets_info[model] = filter_segments_or_datasets(reader.SNP_path, states,
                                                                test_percentiles_list)

            heatmap_data_dir = os.path.join(heatmap_data_path, model + '_tables/')
            if not os.path.isdir(heatmap_data_dir):
                try:
                    os.mkdir(heatmap_data_dir)
                except OSError as exc:
                    if exc.errno != errno.EEXIST:
                        raise
                    pass
            heatmap_data_file = os.path.join(heatmap_data_dir, file_name)

            # print('reading SNP ' + type)
            number_of_datasets, lab, SNP_objects, aligns, segments_number, sum_cov = reader.read_SNPs(method='cover')

            segment_numbers[model] = segments_number
            if cosmic_names.get(cell_line_name):
                corr_to_objects[model] = correlation_with_cosmic(SNP_objects, mode='normal', method='cover',
                                                                 heatmap_data_file=heatmap_data_file,
                                                                 cell_line_name=cell_line_name,
                                                                 cosmic_names=cosmic_names)
            else:
                corr_to_objects[model] = 'NaN'

        for naive_mode in naive_modes:
            if cosmic_names.get(cell_line_name):
                number_of_datasets, lab, SNP_objects, aligns, segments_number, sum_cov = reader.read_SNPs(
                    method=naive_mode)

                corr_to_objects[naive_mode] = correlation_with_cosmic(SNP_objects, mode='normal',
                                                                      cell_line_name=cell_line_name,
                                                                      cosmic_names=cosmic_names)
            else:
                corr_to_objects[naive_mode] = 'NaN'

        # TODO: add 3-5 neighbours naive
        CGH_objects = reader.read_CGH(cgh_names.get(cell_line_name, ''))
        nearest_cgh_objects = find_nearest_probe_to_SNP(SNP_objects, CGH_objects)

        if cosmic_names.get(cell_line_name):
            corr_to_objects_chip = correlation_with_cosmic(CGH_objects, mode='total', cell_line_name=cell_line_name,
                                                           cosmic_names=cosmic_names)
            corr_to_objects_chip_nearest = correlation_with_cosmic(nearest_cgh_objects, mode='total',
                                                                   cell_line_name=cell_line_name,
                                                                   cosmic_names=cosmic_names)
        else:
            corr_to_objects_chip = 'NaN'
            corr_to_objects_chip_nearest = 'NaN'

        out_line = '\t'.join(map(lambda x: '\t'.join(map(str, x)),

                                 [[cell_line_name, lab, aligns, len(SNP_objects), sum_cov, number_of_datasets,
                                   count_cosmic_segments(cell_line_name=cell_line_name, cosmic_names=cosmic_names)]] +

                                 [[segment_numbers[model], corr_to_objects[model]]
                                  for model in map(lambda x: get_name_by_dir(x, naive_modes), snp_dirs)] +

                                 [[corr_to_objects[naive_mode]]
                                  for naive_mode in naive_modes] +

                                 [[corr_to_objects_chip]] +

                                 [[corr_to_objects_chip_nearest]] +

                                 [quality_scores[model]
                                  for model in map(lambda x: get_name_by_dir(x, naive_modes), snp_dirs)] +

                                 [[json.dumps(datasets_info[model])
                                  for model in map(lambda x: get_name_by_dir(x, naive_modes), snp_dirs)]]

                                 )) + '\n'
        out.write(out_line)


if __name__ == '__main__':
    main(sys.argv[1])
