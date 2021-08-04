import os
from .paths_for_components import alignments_path, badmaps_path, tf_dict_path, \
    cl_dict_path, configs_path, results_path, badmaps_dict_path


def get_ending(for_what):
    if for_what == "vcf":
        return ".vcf.gz"
    if for_what == "annotation":
        return ".table_annotated"
    if for_what == "p-value":
        return ".table_p"
    if for_what == "BAD":
        return ".table_BAD"
    if for_what == 'base':
        return ''
    raise AssertionError('Incorrect input parameter')


def create_path_from_master_list_df(row, for_what='base'):
    return os.path.join(alignments_path, row['#EXP'], row['ALIGNS'] + get_ending(for_what))


def get_release_stats_path():
    return os.path.join(results_path, 'release_stats')


def get_result_dir_path(what_for):
    return os.path.join(results_path, '{}_P-values'.format(what_for))


def get_result_table_path(what_for, string):
    return os.path.join(get_result_dir_path(what_for), '{}.tsv'.format(string))


def get_correlation_path():
    return os.path.join(badmaps_path, 'Correlation')


def get_heatmap_data_path():
    return os.path.join(badmaps_path, 'HeatmapData')


def get_badmaps_path_by_validity(valid=False):
    return os.path.join(badmaps_path, 'valid' if valid else 'raw')


def create_badmaps_path_function(name, valid=False):
    return os.path.join(get_badmaps_path_by_validity(valid), 'CAIC', name + ".badmap.tsv")


def create_merged_vcf_path_function(name):
    return os.path.join(badmaps_path, 'merged_vcfs', name + ".tsv")


def get_excluded_badmaps_list_path(remake=False):
    return os.path.join(get_release_stats_path(), 'excluded_badmaps_{}.tsv'.format(2 if remake else 1))


def get_correlation_file_path(remake=False):
    return os.path.join(get_correlation_path(), 'cor_stats{}.tsv'.format('' if remake else '_test'))


def get_new_badmaps_dict_path():
    return os.path.join(os.path.dirname(badmaps_dict_path), 'bad_datasets_dict.json')


def get_merged_badmaps_dict_path(remade=True):
    if remade:
        return os.path.join(os.path.dirname(badmaps_dict_path), 'complete_badmaps_dict.json')
    else:
        return badmaps_dict_path


def get_sarus_dir():
    return os.path.join(results_path, 'Sarus')


def get_sarus_ext(for_what):
    if for_what == 'tsv':
        return '.tsv'
    elif for_what == 'sarus':
        return '.sarus'
    elif for_what == 'fasta':
        return '.fasta'
    else:
        raise AssertionError('Wrong param for get_sarus_ext function {}'.format(for_what))


def get_tf_sarus_path(tf_name, for_what='tsv'):
    return os.path.join(get_sarus_dir(), tf_name + get_sarus_ext(for_what))


def create_neg_bin_stats_path_function(BAD, suffix=''):
    return os.path.join(get_release_stats_path(), 'bias_stats_BAD{:.1f}{}.tsv'.format(
            BAD if BAD else 0, suffix))


def create_neg_bin_weights_path_function(fixed_allele, BAD):
    return os.path.join(configs_path, 'NBweights_{}_BAD={:.1f}.npy'.format(fixed_allele, BAD))


def get_aggregation_dict_path(what_for):
    aggregation_dict_path = None
    if what_for == "TF":
        aggregation_dict_path = tf_dict_path
    if what_for == "CL":
        aggregation_dict_path = cl_dict_path
    if aggregation_dict_path is None:
        raise ValueError("Incorrect usage of open_aggregation_dict function")
    return aggregation_dict_path
