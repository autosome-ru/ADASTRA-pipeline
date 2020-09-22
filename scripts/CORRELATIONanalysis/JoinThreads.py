import os

Correlation_path = '/home/abramov/Correlation_val/'
out_path = Correlation_path + 'cor_stats_test.tsv'


def get_name_by_dir(dir_name):
    return dir_name[:dir_name.rfind('_')].split('/')[-1]


snp_dirs = []
for file_name in sorted(os.listdir(Correlation_path)):
    if file_name.endswith('_tables') and os.path.isdir(Correlation_path + file_name):
        snp_dirs.append(Correlation_path + file_name + '/')

with open(out_path, 'w') as out:
    out.write('\t'.join(map(lambda x: '\t'.join(x),
                            [['#cell_line', 'cells', 'aligns', 'total_snps', 'sum_cov', '#_of_merged_datasets',
                              'total_regions']] +
                            [['number_of_segments_' + get_name_by_dir(snp_dir),
                              'cor_by_snp_' + get_name_by_dir(snp_dir)]
                             for snp_dir in snp_dirs] +
                            [['cor_by_snp_naive',
                              'cor_by_probe_CGH', 'cor_by_snp_probe_CGH']]
                            )) + '\n')
    
    for file_name in os.listdir(Correlation_path):
        if file_name.endswith('.thread') and os.path.isfile(Correlation_path + file_name):
            with open(Correlation_path + file_name) as file:
                out.write(file.readline())
            os.remove(Correlation_path + file_name)
