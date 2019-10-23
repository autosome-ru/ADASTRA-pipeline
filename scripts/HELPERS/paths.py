project_path = "/home/abramov/ASB-Project"
alignments_path = "/home/abramov/Alignments/"

parameters_path = "/home/abramov/PARAMETERS/"
parallel_parameters_path = '/home/abramov/ParallelParameters/'

ploidy_path = "/home/abramov/PloidyForRelease/"
results_path = "/home/abramov/DATA/"

ploidy_dict_path = parameters_path + "CELL_LINES.json"
GTRD_slice_path = parameters_path + "Master-lines.tsv"


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
    return ploidy_path + "Corrected-1,5/" + string + "_ploidy.tsv"
