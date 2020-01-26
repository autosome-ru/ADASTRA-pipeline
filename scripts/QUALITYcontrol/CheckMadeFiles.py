import os.path
import json
import sys
import pandas as pd

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import create_path_from_GTRD_function, make_black_list
from scripts.HELPERS.paths_for_components import results_path, ploidy_dict_path, GTRD_slice_path, tf_dict_path, \
    cl_dict_path, parameters_path


def create_path_for_agr_name(string, agr_name):
    return results_path + agr_name + '_P-values/{}.tsv'.format(string)


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
dict_overall_statistics = {"SNP_calls": None, "unique_SNPs": None, "unique_asb": None, "datasets": None}
for key in dict_overall_statistics:
    dict_overall_statistics[key] = {"TF": {}, "CL": {}}
black_list = make_black_list()
for line in master_list:
    if line[0] == "#":
        continue
    line = line.split("\t")
    if line[0] not in black_list:
        not_blacklisted_exps += 1
        vcf_path = create_path_from_GTRD_function(line, for_what="vcf")
        if os.path.isfile(vcf_path):
            if line[1] not in dict_overall_statistics["datasets"]["TF"]:
                dict_overall_statistics["datasets"]["TF"][line[1]] = 0
            if line[4] not in dict_overall_statistics["datasets"]["CL"]:
                dict_overall_statistics["datasets"]["CL"][line[4]] = 0
            dict_overall_statistics["datasets"]["CL"][line[4]] += 1
            dict_overall_statistics["datasets"]["TF"][line[1]] += 1
            made_experiment_vcfs += 1

        annotated_table_path = create_path_from_GTRD_function(line, for_what="annotated_table")
        if os.path.isfile(annotated_table_path):
            made_annotated_tables += 1

            an_table = pd.read_table(annotated_table_path)
            local_counter = len(an_table.index)
            if local_counter != 0:
                if line[1] not in dict_overall_statistics["SNP_calls"]["TF"]:
                    dict_overall_statistics["SNP_calls"]["TF"][line[1]] = 0
                if line[4] not in dict_overall_statistics["SNP_calls"]["CL"]:
                    dict_overall_statistics["SNP_calls"]["CL"][line[4]] = 0
                dict_overall_statistics["SNP_calls"]["CL"][line[4]] += local_counter
                dict_overall_statistics["SNP_calls"]["TF"][line[1]] += local_counter

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

tf_vcfs_counter = 0
tf_counter = 0
for tf in made_tfs:
    print(tf)
    if os.path.isfile(create_path_for_agr_name(tf, "TF")):
        tf_counter += 1
        for vcf_file in made_tfs[tf]:
            print(vcf_file)
            exp_name = vcf_file.split("/")[-2]
            if os.path.isfile(vcf_file) and exp_name not in black_list:
                tf_vcfs_counter += 1
                tf_table = pd.read_table(vcf_file)

                local_counter = len(tf_table.index)
                fdr_counter = len(tf_table[(tf_table['fdrp_bh_ref'] <= 0.05) | (tf_table["fdrp_bh_alt"] <= 0.05)].index)
                if local_counter != 0:
                    if tf not in dict_overall_statistics["unique_SNPs"]["TF"]:
                        dict_overall_statistics["unique_SNPs"]["TF"][tf] = 0
                    dict_overall_statistics["unique_SNPs"]["TF"][tf] += local_counter

                    if tf not in dict_overall_statistics["unique_asb"]["TF"]:
                        dict_overall_statistics["unique_asb"]["TF"][tf] = 0
                    dict_overall_statistics["unique_asb"]["TF"][tf] += fdr_counter

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

                cl_table = pd.read_table(vcf_file)
                local_counter = len(cl_table.index)
                fdr_counter = len(cl_table[(cl_table['fdrp_bh_ref'] <= 0.05) | (cl_table["fdrp_bh_alt"] <= 0.05)].index)
                if local_counter != 0:
                    if cl not in dict_overall_statistics["unique_SNPs"]["CL"]:
                        dict_overall_statistics["unique_SNPs"]["CL"][cl] = 0
                    dict_overall_statistics["unique_SNPs"]["CL"][cl] += local_counter

                    if cl not in dict_overall_statistics["unique_asb"]["CL"]:
                        dict_overall_statistics["unique_asb"]["CL"][cl] = 0
                    dict_overall_statistics["unique_asb"]["CL"][cl] += fdr_counter


with open(parameters_path + "overall_statistics.JSON", "w") as jsonFile:
    json.dump(dict_overall_statistics, jsonFile)
