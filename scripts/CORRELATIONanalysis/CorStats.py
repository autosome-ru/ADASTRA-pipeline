import os
import sys

sys.path.insert(1, '/home/abramov/ASB-Project')
from scripts.HELPERS.helpers import Reader

scripts_path = '/home/abramov/ASB-Project/scripts/'
Correlation_path = '/home/abramov/Correlation/'

CGH_path = Correlation_path + 'CHIP_hg38.bed'
cosmic_path = Correlation_path + 'COSMIC_copy_number.csv'
synonims_path = scripts_path + 'CORRELATIONanalysis/synonims.tsv'
heatmap_data_path = '/home/abramov/HeatmapData/'


def get_name_by_dir(dir_name):
    if dir_name in naive_names:
        return dir_name
    return dir_name[:dir_name.rfind('_')].split('/')[-1]


if __name__ == '__main__':
    file_name = sys.argv[1]

    out_path = Correlation_path + file_name + '.thread'

    snp_dirs = []
    naive_names = ['naive']

    for f_name in os.listdir(Correlation_path):
        if f_name.endswith('_tables') and os.path.isdir(Correlation_path + f_name):
            snp_dirs.append(Correlation_path + f_name + '/')

    reader = Reader()
    reader.CGH_path = CGH_path
    reader.Cosmic_path = cosmic_path
    reader.synonims_path = synonims_path

    cosmic_names, cgh_names = reader.read_synonims()

    with open(out_path, 'w') as out:

        corr_to_objects_global = dict()
        corr_to_segments_global = dict()

        # if file_name != 'HCT-116_colon_carcinoma_19.tsv': continue

        corr_to_objects = dict()
        corr_to_segments = dict()
        seg_segs = dict()

        # print('reading COSMIC')
        name = file_name[:file_name.rfind('_')]
        index = file_name[file_name.rfind('_') + 1:file_name.rfind('.')]
        COSMIC_segments = reader.read_Cosmic(cosmic_names[name])

        for snp_dir in snp_dirs + naive_names:
            model = get_name_by_dir(snp_dir)
            if model != snp_dir:
                reader.SNP_path = snp_dir + file_name
                method = 'normal'
                cosm_dir = heatmap_data_path + model + '_tables/'
                if not os.path.isdir(cosm_dir):
                    os.mkdir(cosm_dir)
                cosm_path = cosm_dir + file_name
            else:
                method = model
            # print('reading SNP ' + type)
            number_of_datasets, lab, SNP_objects, aligns, segments_number = reader.read_SNPs(method=method)

            model = get_name_by_dir(snp_dir)

            corr_to_objects[model] = 'nan'
            corr_to_segments[model] = 'nan'
            seg_segs[model] = segments_number

            result = COSMIC_segments.merge_with_object_table(SNP_objects)
            if len(result) != 0:
                corr_to_segments[model] = COSMIC_segments.correlation_of_merged(method='mean', result=result)
            if method == 'normal':
                result = SNP_objects.print_merged_to_file(COSMIC_segments, cosm_path)
            else:
                result = SNP_objects.merge_with_segment_table(COSMIC_segments)
            if len(result) != 0:
                corr_to_objects[model] = SNP_objects.correlation_of_merged(result=result)

        # TODO: add 3-5 neighbours naive

        # TODO: add closest chip COR

        # print('reading COSMIC total')
        COSMIC_segments_total = reader.read_Cosmic(cosmic_names[name], mode='total')

        # print('reading CGH')
        N_CGH, CGH_objects = reader.read_CGH(cgh_names[name])

        corr_to_objects_global[name] = 'nan'
        corr_to_segments_global[name] = 'nan'

        result = COSMIC_segments_total.merge_with_object_table(CGH_objects)
        if len(result) != 0:
            corr_to_objects_global[name] = COSMIC_segments_total.correlation_of_merged(method='mean',
                                                                                       result=result)
        result = CGH_objects.merge_with_segment_table(COSMIC_segments_total)
        if len(result) != 0:
            corr_to_segments_global[name] = CGH_objects.correlation_of_merged(result=result)

        out_line = '\t'.join(map(lambda x: '\t'.join(map(str, x)),
                                 [[name, lab, aligns, len(SNP_objects), number_of_datasets, len(COSMIC_segments.segments)]] +
                                 [[seg_segs[model],
                                   corr_to_segments[model],
                                   corr_to_objects[model]]
                                  for model in map(get_name_by_dir, snp_dirs)] +
                                 [[corr_to_segments[name],
                                   corr_to_objects[name]]
                                  for name in naive_names] +
                                 [[corr_to_segments_global[name],
                                   corr_to_objects_global[name]]]
                                 )) + '\n'
        print(file_name)
        out.write(out_line)
