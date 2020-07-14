import sys

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths_for_components import alignments_path, ploidy_path, blacklisted_exps_path, tf_dict_path, \
    cl_dict_path


def create_path_from_GTRD_function(line, for_what, ctrl=False):
    end = ""
    if for_what == "vcf":
        end = ".vcf.gz"
    if for_what == "annotated_table":
        end = "_table_annotated.txt"
    if for_what == "p-value_table":
        end = "_table_p.txt"
    if ctrl:
        return alignments_path + "CTRL/" + line[10] + "/" + line[14] + end
    else:
        return alignments_path + "EXP/" + line[1] + "/" + line[0] + "/" + line[6] + end


def create_ploidy_path_function(string, model):
    return ploidy_path + model + "/" + string + "_ploidy.tsv"


def make_black_list():
    with open(blacklisted_exps_path) as bl:
        black_list = set()
        for line in bl:
            exp_name = line.split(";")[0]
            black_list.add(exp_name)
    return black_list


def open_aggregation_dict(what_for):
    aggregation_dict_path = None
    if what_for == "TF":
        aggregation_dict_path = tf_dict_path
    if what_for == "CL":
        aggregation_dict_path = cl_dict_path
    if aggregation_dict_path is None:
        raise ValueError("Incorrect usage of open_aggregation_dict function")
    return aggregation_dict_path

