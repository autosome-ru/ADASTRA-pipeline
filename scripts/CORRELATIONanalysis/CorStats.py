import os
import sys
import numpy as np
from scipy.stats import kendalltau
from sklearn import linear_model

sys.path.insert(1, '/home/abramov/ASB-Project')
from scripts.HELPERS.helpers import CorrelationReader, Intersection, pack, read_synonims, ChromPos
from scripts.HELPERS.paths import parameters_path, correlation_path, heatmap_data_path

CGH_path = parameters_path + 'CHIP_hg38.sorted.bed'
cosmic_path = parameters_path + 'COSMIC_copy_number.sorted.tsv'


def get_name_by_dir(dir_name):
    if dir_name in naive_modes:
        return dir_name
    return dir_name[:dir_name.rfind('_')].split('/')[-1]


def count_cosmic_segments():
    with open(cosmic_path, 'r') as cosmic_file:
        return sum(True for line in cosmic_file if unpack_cosmic_segments(line))


def unpack_cosmic_segments(line, mode='normal'):
    if line[0] == '#':
        return [''] * len(line.strip().split('\t'))
    line = line.strip().split('\t')

    if line[0] != cosmic_names[cell_line_name]:
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


def correlation_with_cosmic(SNP_objects, mode, heatmap_data_file=None):
    if heatmap_data_file is not None:
        heatmap = open(heatmap_data_file, 'w')
    cosmic_segments = []
    with open(cosmic_path, 'r') as cosmic_file:
        for line in cosmic_file:
            # TODO: change v name
            v = unpack_cosmic_segments(line, mode=mode)
            if v:
                cosmic_segments.append(v)
    snp_ploidy = []
    cosm_ploidy = []
    for chr, pos, ploidy, qual, segn, in_intersect, cosmic_ploidy \
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
        lm = linear_model.LinearRegression()
        lm.fit(np.array(snp_ploidy), np.array(cosm_ploidy))
        return kendalltau(snp_ploidy, cosm_ploidy)[0], tuple(map(float, lm.predict(np.array([0, 1]).reshape(-1, 1))))
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


if __name__ == '__main__':
    file_name = sys.argv[1]
    print(file_name)

    out_path = correlation_path + file_name + '.thread'

    snp_dirs = []
    naive_modes = ['naive']

    for f_name in os.listdir(correlation_path):
        if f_name.endswith('_tables') and os.path.isdir(correlation_path + f_name):
            snp_dirs.append(correlation_path + f_name + '/')

    reader = CorrelationReader()
    reader.CGH_path = CGH_path

    cosmic_names, cgh_names = read_synonims()

    with open(out_path, 'w') as out:

        # if file_name != 'HCT-116_colon_carcinoma_19.tsv': continue

        corr_to_objects = {}
        segment_numbers = {}
        lm_coefficients = {}
        # print('reading COSMIC')
        cell_line_name = file_name[:file_name.rfind('_')]
        index = file_name[file_name.rfind('_') + 1:file_name.rfind('.')]

        for snp_dir in snp_dirs:
            model = get_name_by_dir(snp_dir)

            reader.SNP_path = snp_dir + file_name

            heatmap_data_dir = heatmap_data_path + model + '_tables/'
            if not os.path.isdir(heatmap_data_dir):
                os.mkdir(heatmap_data_dir)
            heatmap_data_file = heatmap_data_dir + file_name

            # print('reading SNP ' + type)
            number_of_datasets, lab, SNP_objects, aligns, segments_number = reader.read_SNPs(method='normal')

            segment_numbers[model] = segments_number
            corr_to_objects[model], lm_coefficients[model] = correlation_with_cosmic(
                                                                                    SNP_objects,
                                                                                    mode='normal',
                                                                                    heatmap_data_file=heatmap_data_file
            )

        for naive_mode in naive_modes:
            number_of_datasets, lab, SNP_objects, aligns, segments_number = reader.read_SNPs(method=naive_mode)

            corr_to_objects[naive_mode], _ = correlation_with_cosmic(SNP_objects, mode='normal')

        # TODO: add 3-5 neighbours naive
        CGH_objects = reader.read_CGH(cgh_names[cell_line_name])
        nearest_cgh_objects = find_nearest_probe_to_SNP(SNP_objects, CGH_objects)

        corr_to_objects_chip, _ = correlation_with_cosmic(CGH_objects, mode='total')
        corr_to_objects_chip_nearest, _ = correlation_with_cosmic(nearest_cgh_objects, mode='total')
        out_line = '\t'.join(map(lambda x: '\t'.join(map(str, x)),

                                 [[cell_line_name, lab, aligns, len(SNP_objects), number_of_datasets,
                                   count_cosmic_segments()]] +

                                 [[segment_numbers[model], corr_to_objects[model],
                                   lm_coefficients[model][0],lm_coefficients[model[1]]]
                                  for model in map(get_name_by_dir, snp_dirs)] +

                                 [[corr_to_objects[naive_mode]]
                                  for naive_mode in naive_modes] +

                                 [[corr_to_objects_chip]] +

                                 [[corr_to_objects_chip_nearest]]

                                 )) + '\n'
        out.write(out_line)
