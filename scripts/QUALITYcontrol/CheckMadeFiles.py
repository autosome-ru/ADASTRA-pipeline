import os.path
import json
import sys

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import ploidy_dict_path, create_path_from_GTRD_function, results_path, \
    GTRD_slice_path, create_ploidy_path_function, make_black_list, cl_dict_path, tf_dict_path


def create_path_for_agr_name(string, agr_name):
    return results_path + agr_name + "_P-values/" + string + '_common_table.tsv'


with open(ploidy_dict_path, "r") as cl_file, \
        open(GTRD_slice_path, "r") as ml, \
        open(tf_dict_path, "r") as tfs, \
        open(cl_dict_path, "r") as cls:
    cell_lines = json.loads(cl_file.readline())
    made_tfs = json.loads(tfs.readline())
    made_cls = json.loads(cls.readline())
    master_list = ml.readlines()

not_blacklisted_exps = 0
not_blacklisted_ctrl = 0
made_experiment_vcfs = 0
made_control_vcfs = 0
made_p_tables = 0
made_annotated_tables = 0
counted_controls = set()
SNP_counter = 0
dict_SNP_TF_statistics = {}

black_list = make_black_list()
for line in master_list:
    if line[0] == "#":
        continue
    line = line.split("\t")
    if line[0] not in black_list:
        not_blacklisted_exps += 1
        vcf_path = create_path_from_GTRD_function(line, for_what="vcf")
        if os.path.isfile(vcf_path):
            made_experiment_vcfs += 1

        p_value_table_path = create_path_from_GTRD_function(line, for_what="p-value_table")
        if os.path.isfile(p_value_table_path):
            made_p_tables += 1

        annotated_table_path = create_path_from_GTRD_function(line, for_what="annotated_table")
        if os.path.isfile(annotated_table_path):
            made_annotated_tables += 1
            with open(annotated_table_path, "r") as an_table:
                local_counter = 0
                for SNP in an_table:
                    if SNP[0] == "#":
                        continue
                    local_counter += 1
                    SNP_counter += 1
                if local_counter != 0:
                    try:
                        el = dict_SNP_TF_statistics[line[1]]
                        dict_SNP_TF_statistics[line[1]] = el + local_counter
                    except KeyError:
                        dict_SNP_TF_statistics[line[1]] = local_counter

        if len(line) > 10 and line[10] not in black_list:
            vcf_path = create_path_from_GTRD_function(line, for_what="vcf", ctrl=True)
            if vcf_path in counted_controls:
                continue
            not_blacklisted_ctrl += 1
            counted_controls.add(vcf_path)
            if os.path.isfile(vcf_path):
                made_control_vcfs += 1

print("Made {}/{} VCFS ({}/{} experiment VCFs, {}/{} control VCFs), {} annotated tables, {} P-value tables".format(
    made_control_vcfs + made_experiment_vcfs, not_blacklisted_ctrl + not_blacklisted_exps,
    made_experiment_vcfs, not_blacklisted_exps,
    made_control_vcfs, not_blacklisted_ctrl,
    made_annotated_tables, made_p_tables))

print("Total of {} SNPs in experiment VCFs".format(SNP_counter))
SNP_TF_statistics = sorted(list(dict_SNP_TF_statistics.items()), key=lambda x: x[1], reverse=True)
print(SNP_TF_statistics[:5])

ploidy_control_vcfs = 0
ploidy_vcfs_counter = 0
ploidy_counter = 0
counted_control_vcfs = set()
for ploidy in cell_lines:
    ploidy_file = create_ploidy_path_function(ploidy, model='CAIC')
    if os.path.isfile(ploidy_file):
        ploidy_counter += 1
        for vcf_file in cell_lines[ploidy]:
            exp_name = vcf_file.split("/")[-2]
            if os.path.isfile(vcf_file) and exp_name not in black_list:
                ploidy_vcfs_counter += 1
                if vcf_file.find("CTRL") != -1 and vcf_file not in counted_control_vcfs:
                    ploidy_control_vcfs += 1
                    counted_control_vcfs.add(vcf_file)
print("Made {} ploidies from  {} VCFs ({} experiment VCFs)".format(ploidy_counter, ploidy_vcfs_counter,
                                                                   ploidy_vcfs_counter - ploidy_control_vcfs))

tf_vcfs_counter = 0
tf_counter = 0
for tf in made_tfs:
    if os.path.isfile(create_path_for_agr_name(tf, "TF")):
        tf_counter += 1
        for vcf_file in made_tfs[tf]:
            exp_name = vcf_file.split("/")[-2]
            if os.path.isfile(vcf_file) and exp_name not in black_list:
                tf_vcfs_counter += 1
print("Made aggregation for {} TFs from  {} VCFs".format(tf_counter, tf_vcfs_counter))


cl_vcfs_counter = 0
cl_counter = 0
for cl in made_cls:
    if os.path.isfile(create_path_for_agr_name(cl, "CL")):
        cl_counter += 1
        for vcf_file in made_cls[cl]:
            exp_name = vcf_file.split("/")[-2]
            if os.path.isfile(vcf_file) and exp_name not in black_list:
                cl_vcfs_counter += 1
print("Made aggregation for {} cell lines from  {} VCFs".format(cl_counter, cl_vcfs_counter))
