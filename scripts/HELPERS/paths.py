import os

from .paths_for_components import alignments_path, badmaps_path, tf_dict_path, \
    cl_dict_path, configs_path, results_path, badmaps_dict_path

stage_dict = {
    'BAD': 'BAD_annotations',
    'p-value': 'Pvalue_annotations',
    'Sarus': 'Sarus'
}


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
    if for_what in stage_dict:
        return os.path.join(get_dir_by_stage(for_what), row['ALIGNS'] + get_ending(for_what))
    else:
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


def create_merged_vcf_path_function(name):
    return os.path.join(badmaps_path, 'merged_vcfs', name + ".tsv")


def get_excluded_badmaps_list_path(model, remake=False):  #FIXME
    return os.path.join(get_release_stats_path(), 'excluded_badmaps_{}_{}.tsv'.format(2 if remake else 1, model))


def get_correlation_file_path(remake=False):
    return os.path.join(get_correlation_path(), 'cor_stats{}.tsv'.format('' if remake else '_test'))


def get_new_badmaps_dict_path(model):
    return os.path.join(os.path.dirname(badmaps_dict_path), 'bad_datasets_dict_{}.json'.format(model))


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
    bad_dir = os.path.join(get_release_stats_path(),
                           'BAD{:.2f}'.format(BAD if BAD else 0))
    if not os.path.exists(bad_dir):
        os.mkdir(bad_dir)
    return os.path.join(bad_dir, 'stats{}.tsv'.format(suffix))


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


def get_dir_by_stage(stage):
    return os.path.join(results_path, stage_dict[stage])
