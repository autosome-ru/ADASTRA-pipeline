import os

from scripts.HELPERS.helpers import test_percentiles_list
from scripts.HELPERS.paths import get_correlation_path


correlation_path = get_correlation_path()
out_path = os.path.join(correlation_path, 'cor_stats_test.tsv')


def get_name_by_dir(dir_name):
    return dir_name[:dir_name.rfind('_')].split('/')[-1]


def main():
    snp_dirs = []
    for file_name in sorted(os.listdir(correlation_path)):
        snp_dir = os.path.join(correlation_path, file_name)
        if file_name.endswith('_tables') and os.path.isdir(snp_dir):
            snp_dirs.append(snp_dir)

    with open(out_path, 'w') as out:
        out.write('\t'.join(map(lambda x: '\t'.join(x),
                                [['#cell_line', 'cells', 'aligns', 'total_snps', 'sum_cov', '#_of_merged_datasets',
                                  'total_regions']] +
                                [['number_of_segments_' + get_name_by_dir(snp_dir),
                                  'cor_by_snp_' + get_name_by_dir(snp_dir)]
                                 for snp_dir in snp_dirs] +
                                [['cor_by_snp_naive',
                                  'cor_by_probe_CGH', 'cor_by_snp_probe_CGH']] +
                                [['Q{}_{}'.format(percentile, get_name_by_dir(snp_dir)) for
                                  percentile in test_percentiles_list]
                                 for snp_dir in snp_dirs]
                                )) + '\n')

        for file_name in os.listdir(correlation_path):
            thread = os.path.join(correlation_path, file_name)
            if file_name.endswith('.thread') and os.path.isfile(thread):
                with open(thread) as file:
                    out.write(file.readline())
                os.remove(thread)


if __name__ == '__main__':
    main()
