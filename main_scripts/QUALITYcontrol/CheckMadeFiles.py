import os.path
import json

parameters_path = "/home/abramov/PARAMETERS/"
alignments_path = "/home/abramov/Alignments/"


def create_path(line, for_what, ctrl=False):
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


def make_black_list():
    with open(parameters_path + "blacklisted_exps.tsv") as bl:
        black_list = set()
        for line in bl:
            exp_name = line.split("\t")[0]
            black_list.add(exp_name)
    return black_list


with open(parameters_path + "CELL_LINES.json", "r") as f, open(parameters_path + "Master-lines.tsv", "r") as ml:
    cell_lines = json.loads(f.readline())
    master_list = ml.readlines()

made_experiment_vcfs = 0
made_control_vcfs = 0
made_p_tables = 0
made_annotated_tables = 0

black_list = make_black_list()
for line in master_list:
    line = line.split("\t")
    if line[0] not in black_list:

        vcf_path = create_path(line, for_what="vcf")
        if os.path.isfile(vcf_path):
            made_experiment_vcfs += 1

        p_value_table_path = create_path(line, for_what="p-value_table")
        if os.path.isfile(p_value_table_path):
            made_p_tables += 1

        annotated_table_path = create_path(line, for_what="annotated_table")
        if os.path.isfile(annotated_table_path):
            made_annotated_tables += 1

        if len(line) > 10 and line[10] not in black_list:
            vcf_path = create_path(line, for_what="vcf", ctrl=True)
            if os.path.isfile(vcf_path):
                made_control_vcfs += 1
print("Made {} experiment VCFs, {} control VCFs, {} annotated tables, {} P-value tables".format(made_experiment_vcfs,
                                                                                                made_control_vcfs,
                                                                                                made_annotated_tables,
                                                                                                made_p_tables,))
