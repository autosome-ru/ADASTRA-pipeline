import os.path
import json

parameters_path = "/home/abramov/PARAMETERS/"
alignments_path = "/home/abramov/Alignments/"
ploidy_path = "/home/abramov/PloidyForRelease/"

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


def create_ploidy(string):
    path = ploidy_path + "Corrected-1,5/" + string + "_ploidy.tsv"
    return path


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
counted_controls = set()

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
            if os.path.isfile(vcf_path) and vcf_path not in counted_controls:
                made_control_vcfs += 1
                counted_controls.add(vcf_path)
print("Made {} experiment VCFs, {} control VCFs, {} annotated tables, {} P-value tables".format(made_experiment_vcfs,
                                                                                                made_control_vcfs,
                                                                                                made_annotated_tables,
                                                                                                made_p_tables,))
ploidy_control_vcfs = 0
ploidy_vcfs_counter = 0
ploidy_counter = 0
counted_control_vcfs = set()
for ploidy in cell_lines:
    print(ploidy)
    ploidy_file = create_ploidy(ploidy)
    if os.path.isfile(ploidy_file):
        ploidy_counter += 1
        for vcf_file in cell_lines[ploidy]:
            exp_name = vcf_file.split("/")[-2:-1]
            if os.path.isfile(vcf_file) and exp_name not in black_list:
                ploidy_vcfs_counter += 1
                if vcf_file.find("CTRL") != -1 and vcf_file not in counted_control_vcfs:
                    ploidy_control_vcfs += 1
                    counted_control_vcfs.add(vcf_file)
print("Made {} ploidies from  {} VCFs ({} experiment VCFs)".format(ploidy_counter, ploidy_vcfs_counter,
                                                                   ploidy_vcfs_counter - ploidy_control_vcfs))
