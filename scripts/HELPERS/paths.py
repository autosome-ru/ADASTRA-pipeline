import os
from scripts.HELPERS.paths_for_components import alignments_path, badmaps_path, tf_dict_path, \
    cl_dict_path


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


def create_badmaps_path_function(name):
    return os.path.join(badmaps_path, 'CAIC', name + ".badmap.tsv")


def create_merged_vcf_path_function(name):
    return os.path.join(badmaps_path, 'merged_vcfs', name + ".tsv")


def open_aggregation_dict(what_for):
    aggregation_dict_path = None
    if what_for == "TF":
        aggregation_dict_path = tf_dict_path
    if what_for == "CL":
        aggregation_dict_path = cl_dict_path
    if aggregation_dict_path is None:
        raise ValueError("Incorrect usage of open_aggregation_dict function")
    return aggregation_dict_path
