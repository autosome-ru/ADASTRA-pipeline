import os
import sys
from scipy.stats import kendalltau

sys.path.insert(1, '/home/abramov/ASB-Project')
from scripts.HELPERS.helpers import CorrelationReader, ChromPos, Intersection, pack, read_synonims
from scripts.HELPERS.paths import parameters_path, correlation_path, heatmap_data_path

CGH_path = parameters_path + 'CHIP_hg38.bed'
cosmic_path = parameters_path + 'COSMIC_copy_number.csv'


def get_name_by_dir(dir_name):
    if dir_name in naive_modes:
        return dir_name
    return dir_name[:dir_name.rfind('_')].split('/')[-1]


def count_cosmic_segments():
    with open(cosmic_path, 'r') as cosmic_file:
        return sum(True for line in cosmic_file if unpack_cosmic_segments(line))


def unpack_cosmic_segments(line, mode='normal'):
    if line[0] == '#':
        return []
    line = line.strip().split(',')
    if line[0] != name:
        return []
    if 'chr' + line[4] not in ChromPos.chrs:
        return []
    if int(line[10]) == 0:
        return []
    
    if mode == 'normal':
        value = int(line[11]) / int(line[10]) - 1
    elif mode == 'total':
        value = int(line[11])
    else:
        raise ValueError(mode)
    
    return ['chr' + line[4], int(line[5]), int(line[6]), value]


def correlation_with_cosmic(SNPs_iterator, mode, heatmap_file_path=None):
    heat_map = None
    if heatmap_file_path is not None:
        heat_map = open(heatmap_file_path, 'w')
    with open(cosmic_path, 'r') as cosmic_file:
        snp_bad_list = []
        cosmic_bad_list = []
        for chr, pos, ploidy, qual, segn, in_intersect, cosmic_bad \
                in Intersection(SNPs_iterator, cosmic_file, write_intersect=True,
                                unpack_segments_function=lambda x: unpack_cosmic_segments(x, mode=mode),
                                write_segment_args=True):
            if not in_intersect:
                continue
            snp_bad_list.append(ploidy)
            cosmic_bad_list.append(cosmic_bad)

            if heat_map is not None:
                heat_map.write(pack([chr, pos, ploidy, cosmic_bad]))
    if heat_map is not None:
        heat_map.close()
    
    if len(snp_bad_list) != 0:
        return kendalltau(snp_bad_list, cosmic_bad_list)[0]
    return 'NaN'


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
        
        corr_to_objects = dict()
        segment_numbers = dict()
        
        # print('reading COSMIC')
        name = file_name[:file_name.rfind('_')]
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
            corr_to_objects[model] = correlation_with_cosmic(SNP_objects,
                                                             mode='normal',
                                                             heatmap_file_path=heatmap_data_file)
        
        for naive_mode in naive_modes:
            number_of_datasets, lab, SNP_objects, aligns, segments_number = reader.read_SNPs(method=naive_mode)
            
            corr_to_objects[model] = correlation_with_cosmic(SNP_objects, mode='normal')
        
        # TODO: add 3-5 neighbours naive
        
        # TODO: add closest chip COR
        
        # print('reading CGH')
        CGH_objects = reader.read_CGH(cgh_names[name])
        
        corr_to_objects_chip = correlation_with_cosmic(CGH_objects, mode='total')
        
        out_line = '\t'.join(map(lambda x: '\t'.join(map(str, x)),
        
                                 [[name, lab, aligns, len(SNP_objects), number_of_datasets,
                                   count_cosmic_segments()]] +
        
                                 [[segment_numbers[model], corr_to_objects[model]]
                                  for model in map(get_name_by_dir, snp_dirs)] +
        
                                 [[corr_to_objects[naive_mode]]
                                  for naive_mode in naive_modes] +
        
                                 [[corr_to_objects_chip]]
        
                                 )) + '\n'
        out.write(out_line)
