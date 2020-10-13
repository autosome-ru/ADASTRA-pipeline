import os.path
import json
import pandas as pd

from scripts.HELPERS.helpers import remove_punctuation, dtype_dict
from scripts.HELPERS.paths import create_path_from_master_list_df, get_result_table_path, get_result_dir_path
from scripts.HELPERS.paths_for_components import results_path, parameters_path, badmaps_dict_path, master_list_path


def main():
    made_experiment_vcfs = 0
    made_control_vcfs = 0

    made_annotated_tables = 0
    dict_overall_statistics = {"SNP_calls": None, "unique_SNPs": None, "unique_asb": None, "datasets": None,
                               "unique_SNPs_rs": None, "unique_asb_rs": None}
    for key in dict_overall_statistics:
        dict_overall_statistics[key] = {"TF": {}, "CL": {}}
    ml_df = pd.read_table(master_list_path, dtype=dtype_dict)
    ml_df['vcf_path'] = ml_df.apply(create_path_from_master_list_df, axis=1)
    ml_df = ml_df[ml_df['vcf_path'].apply(os.path.isfile)]

    for index, row in ml_df.iterrows():
        row['CELLS'] = remove_punctuation(row['CELLS'])
        if row['TF_UNIPROT_ID'] not in dict_overall_statistics["datasets"]["TF"]:
            dict_overall_statistics["datasets"]["TF"][row['TF_UNIPROT_ID']] = 0
        if row['CELLS'] not in dict_overall_statistics["datasets"]["CL"]:
            dict_overall_statistics["datasets"]["CL"][row['CELLS']] = 0
        dict_overall_statistics["datasets"]["CL"][row['CELLS']] += 1
        dict_overall_statistics["datasets"]["TF"][row['TF_UNIPROT_ID']] += 1
        made_experiment_vcfs += 1
        annotated_table_path = create_path_from_master_list_df(row, for_what="annotation")
        if os.path.isfile(annotated_table_path):
            made_annotated_tables += 1
        an_table = pd.read_table(annotated_table_path)
        if an_table.empty:
            continue
        an_table = an_table[(an_table['ID'] != '.') & (an_table['ID'])]
        local_counter = len(an_table.index)
        if local_counter == 0:
            continue
        if row['TF_UNIPROT_ID'] not in dict_overall_statistics["SNP_calls"]["TF"]:
            dict_overall_statistics["SNP_calls"]["TF"][row['TF_UNIPROT_ID']] = 0
        if row['CELLS'] not in dict_overall_statistics["SNP_calls"]["CL"]:
            dict_overall_statistics["SNP_calls"]["CL"][row['CELLS']] = 0
        dict_overall_statistics["SNP_calls"]["CL"][row['CELLS']] += local_counter
        dict_overall_statistics["SNP_calls"]["TF"][row['TF_UNIPROT_ID']] += local_counter

    print("Made {} VCFS ({} experiment VCFs, {} control VCFs), {} annotated tables".format(
        made_control_vcfs + made_experiment_vcfs,
        made_experiment_vcfs,
        made_control_vcfs,
        made_annotated_tables))

    vcf_counter = {'TF': 0, 'CL': 0}
    total_fdrs = 0
    total_fdrs_rs = 0
    for what_for in ('TF', 'CL'):
        for obj in os.listdir(get_result_dir_path(what_for)):
            vcf_counter[what_for] += 1
            obj_table = pd.read_table(get_result_table_path(obj, what_for))
            if obj_table.empty:
                continue
            local_counter = len(obj_table.index)
            local_counter_rs = len(obj_table['ID'].unique())
            fdr_counter = len(obj_table[(obj_table['fdrp_bh_ref'] <= 0.05) | (obj_table["fdrp_bh_alt"] <= 0.05)].index)
            fdr_counter_rs = len(obj_table[(obj_table['fdrp_bh_ref'] <= 0.05) | (obj_table["fdrp_bh_alt"] <= 0.05)]['ID'].unique())
            if obj not in dict_overall_statistics["unique_SNPs"][what_for]:
                dict_overall_statistics["unique_SNPs"][what_for][obj] = 0
                dict_overall_statistics["unique_SNPs_rs"][what_for][obj] = 0
            dict_overall_statistics["unique_SNPs"][what_for][obj] += local_counter
            dict_overall_statistics["unique_SNPs_rs"][what_for][obj] += local_counter_rs

            if obj not in dict_overall_statistics["unique_asb"][what_for]:
                dict_overall_statistics["unique_asb"][what_for][obj] = 0
                dict_overall_statistics["unique_asb_rs"][what_for][obj] = 0
            dict_overall_statistics["unique_asb"][what_for][obj] += fdr_counter
            dict_overall_statistics["unique_asb_rs"][what_for][obj] += fdr_counter_rs
            total_fdrs += fdr_counter
            total_fdrs_rs += fdr_counter_rs
        print("Made aggregation for {} {}s".format(vcf_counter[what_for], what_for))
        print('In {} aggregation total of {} ASB events'.format(what_for, total_fdrs))

    with open(os.path.join(parameters_path, "overall_statistics.json"), "w") as json_file:
        json.dump(dict_overall_statistics, json_file)


if __name__ == '__main__':
    main()
