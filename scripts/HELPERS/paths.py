import sys
sys.path.insert(1, "/home/abramov/ASB-Project")

project_path = "/home/abramov/ASB-Project/"
alignments_path = "/home/abramov/Alignments/"

scripts_path = project_path + "scripts/"
parameters_path = "/home/abramov/PARAMETERS/"
parallel_parameters_path = '/home/abramov/ParallelParameters/stats'

ploidy_path = "/home/abramov/StatsPloidy/"
results_path = "/home/abramov/DATA/"

ploidy_dict_path = parameters_path + "CELL_LINES.json"
GTRD_slice_path = parameters_path + "Master-lines.tsv"
blacklisted_exps_path = parameters_path + "blacklisted_exps.tsv"
tf_dict_path = parameters_path + "TF_DICT.json"
cl_dict_path = parameters_path + "CL_DICT.json"

correlation_path = '/home/abramov/Correlation/'
synonims_path = parameters_path + 'synonims.tsv'
heatmap_data_path = '/home/abramov/HeatmapData/'


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


def create_ploidy_path_function(string):
    return ploidy_path + "Corrected-6/" + string + "_ploidy.tsv"


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

