import os
import re
import sys
from scipy.stats import kendalltau
import numpy as np
import pandas as pd

from scripts.HELPERS.helpers import CorrelationReader, Intersection, pack, read_synonims, ChromPos, get_states
from scripts.HELPERS.paths_for_components import correlation_path, heatmap_data_path, cgh_path, cosmic_path


def get_name_by_dir(dir_name, naive_modes):
    if dir_name in naive_modes:
        return dir_name
    return dir_name[:dir_name.rfind('_')].split('/')[-1]


def count_cosmic_segments(cosmic_names, cell_line_name):
    with open(cosmic_path, 'r') as cosmic_file:
        return sum(True for line in cosmic_file if unpack_cosmic_segments(line, cell_line_name=cell_line_name, cosmic_names=cosmic_names))


def unpack_cosmic_segments(line, mode='normal', cosmic_names='', cell_line_name=''):
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


def correlation_with_cosmic(SNP_objects, mode, method='normal', heatmap_data_file=None, cell_line_name='', cosmic_names=''):
    if heatmap_data_file is not None:
        heatmap = open(heatmap_data_file, 'w')
    cosmic_segments = []
    with open(cosmic_path, 'r') as cosmic_file:
        for line in cosmic_file:
            # TODO: change v name
            v = unpack_cosmic_segments(line, mode=mode, cell_line_name=cell_line_name, cosmic_names=cosmic_names)
            if v:
                cosmic_segments.append(v)
    snp_ploidy = []
    cosm_ploidy = []
    if method == 'normal':
        for chr, pos, ploidy, quals, in_intersect, cosmic_ploidy \
                in Intersection(SNP_objects, cosmic_segments, write_intersect=True,
                                write_segment_args=True):
            if not in_intersect:
                continue
            snp_ploidy.append(ploidy)
            cosm_ploidy.append(cosmic_ploidy)

            if heatmap_data_file is not None:
                heatmap.write(pack([chr, pos, ploidy, cosmic_ploidy]))
        if heatmap_data_file is not None:
            heatmap.close()

        if len(snp_ploidy) != 0:
            kt = kendalltau(snp_ploidy, cosm_ploidy)[0]
            if kt == 'nan':
                return 'NaN'
            return kt
        return 'NaN'
    elif method == 'cover':

        for chr, pos, cov, ploidy, quals, in_intersect, cosmic_ploidy \
                in Intersection(SNP_objects, cosmic_segments, write_intersect=True,
                                write_segment_args=True):
            if not in_intersect:
                continue
            snp_ploidy.append(ploidy)
            cosm_ploidy.append(cosmic_ploidy)

            if heatmap_data_file is not None:
                heatmap.write(pack([chr, pos, cov, ploidy, cosmic_ploidy] + [quals[x] for x in quals]))
        if heatmap_data_file is not None:
            heatmap.close()

        if len(snp_ploidy) != 0:
            kt = kendalltau(snp_ploidy, cosm_ploidy)[0]
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


def filter_segments_or_datasets(snps_path, states):
    with open(snps_path, 'r') as out:
        header_comment = out.readline()
        if not out.readline():
            return
    out_table = pd.read_table(snps_path, header=None, comment='#')
    out_table.columns = ['chr', 'pos', 'ref', 'alt', 'BAD'] + ['Q{:.2f}'.format(BAD) for BAD in states] + ['snps_n',
                                                                                                           'sumcov',
                                                                                                           'dataset',
                                                                                                           'seg_id',
                                                                                                           'p_value']
    valid_segments = set(dataset for dataset in list(set(out_table['dataset'])) if
                         np.quantile(out_table[out_table['dataset'] == dataset]['p_value'], 0.05) >= 0.05 and len(
                             out_table[out_table['dataset'] == dataset].index) >= 10)
    out_table = out_table[out_table['seg_id'].isin(valid_segments)]
    with open(snps_path, 'w') as out:
        out.write(header_comment)
    out_table.to_csv(snps_path, header=False, index=False, sep='\t', mode='a')


def main(file_name):
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
        # print('reading COSMIC')
        cell_line_name = file_name[:file_name.rfind('_')]
        index = file_name[file_name.rfind('_') + 1:file_name.rfind('.')]

        for snp_dir in snp_dirs:
            model = get_name_by_dir(snp_dir, naive_modes)

            if re.match(r'^CAIC@.+@.+$', model) is not None:
                print(model.split('@')[1])
                states = get_states(model.split('@')[1])
            else:
                states = get_states('')
            reader.states = states

            reader.SNP_path = os.path.join(snp_dir, file_name)

            filter_segments_or_datasets(reader.SNP_path, states)

            heatmap_data_dir = os.path.join(heatmap_data_path, model + '_tables/')
            if not os.path.isdir(heatmap_data_dir):
                try:
                    os.mkdir(heatmap_data_dir)
                except:
                    pass
            heatmap_data_file = os.path.join(heatmap_data_dir, file_name)

            # print('reading SNP ' + type)
            number_of_datasets, lab, SNP_objects, aligns, segments_number, sum_cov = reader.read_SNPs(method='cover')

            segment_numbers[model] = segments_number
            if cosmic_names.get(cell_line_name):
                corr_to_objects[model] = correlation_with_cosmic(SNP_objects, mode='normal', method='cover',
                                                                 heatmap_data_file=heatmap_data_file, cell_line_name=cell_line_name, cosmic_names=cosmic_names)
            else:
                corr_to_objects[model] = 'NaN'

        for naive_mode in naive_modes:
            if cosmic_names.get(cell_line_name):
                number_of_datasets, lab, SNP_objects, aligns, segments_number, sum_cov = reader.read_SNPs(
                    method=naive_mode)

                corr_to_objects[naive_mode] = correlation_with_cosmic(SNP_objects, mode='normal', cell_line_name=cell_line_name, cosmic_names=cosmic_names)
            else:
                corr_to_objects[naive_mode] = 'NaN'

        # TODO: add 3-5 neighbours naive
        CGH_objects = reader.read_CGH(cgh_names.get(cell_line_name, ''))
        nearest_cgh_objects = find_nearest_probe_to_SNP(SNP_objects, CGH_objects)

        if cosmic_names.get(cell_line_name):
            corr_to_objects_chip = correlation_with_cosmic(CGH_objects, mode='total', cell_line_name=cell_line_name, cosmic_names=cosmic_names)
            corr_to_objects_chip_nearest = correlation_with_cosmic(nearest_cgh_objects, mode='total', cell_line_name=cell_line_name, cosmic_names=cosmic_names)
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

                                 [[corr_to_objects_chip_nearest]]

                                 )) + '\n'
        out.write(out_line)


if __name__ == '__main__':
    main(sys.argv[1])
