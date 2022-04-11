import os.path
import json
import pandas as pd

from scripts.HELPERS.helpers import remove_punctuation, dtype_dict
from scripts.HELPERS.paths import create_path_from_master_list_df, get_result_table_path, get_result_dir_path, \
    get_release_stats_path
from scripts.HELPERS.paths_for_components import master_list_path


def main():
    made_experiment_vcfs = 0
    made_control_vcfs = 0

    made_annotated_tables = 0
    dict_overall_statistics = {"SNP_calls": None, 'raw_SNP_calls': None,
                               "unique_SNPs": None,
                               "unique_asb": None, "datasets": None,
                               "unique_SNPs_rs": None, "unique_asb_rs": None,
                               }
    for key in dict_overall_statistics:
        dict_overall_statistics[key] = {"TF": {}, "CL": {}}
    ml_df = pd.read_table(master_list_path, dtype=dtype_dict)
    total_vcf_count = len(ml_df.index)
    ml_df['vcf_path'] = ml_df.apply(lambda x: create_path_from_master_list_df(x, 'vcf'), axis=1)
    not_downloaded_df = ml_df[~ml_df['vcf_path'].apply(os.path.isfile)]
    not_downloaded_df.loc[:, not_downloaded_df.columns != 'vcf_path'].to_csv(
        os.path.join(get_release_stats_path(), 'not_downloaded.tsv'),
        index=False, sep='\t')
    ml_df = ml_df[ml_df['vcf_path'].apply(os.path.isfile)]
    for index, row in ml_df.iterrows():
        row['CELLS'] = remove_punctuation(row['CELLS'])
        if row['CELLS'] not in dict_overall_statistics["datasets"]["CL"]:
            dict_overall_statistics["datasets"]["CL"][row['CELLS']] = 0
        dict_overall_statistics["datasets"]["CL"][row['CELLS']] += 1
        made_experiment_vcfs += 1
        vcf_df = pd.read_table(row['vcf_path'], header=None, comment='#')
        if vcf_df.empty:
            continue
        local_counter = len(vcf_df.index)
        if row['CELLS'] not in dict_overall_statistics["raw_SNP_calls"]["CL"]:
            dict_overall_statistics["raw_SNP_calls"]["CL"][row['CELLS']] = 0
        dict_overall_statistics["raw_SNP_calls"]["CL"][row['CELLS']] += local_counter
        annotated_table_path = create_path_from_master_list_df(row, for_what="annotation")
        if not os.path.isfile(annotated_table_path):
            continue
        made_annotated_tables += 1
        an_table = pd.read_table(annotated_table_path)
        if an_table.empty:
            continue
        an_table = an_table[(an_table['ID'] != '.') & (an_table['ID'])]
        local_counter = len(an_table.index)
        if local_counter == 0:
            continue
        if row['CELLS'] not in dict_overall_statistics["SNP_calls"]["CL"]:
            dict_overall_statistics["SNP_calls"]["CL"][row['CELLS']] = 0
        dict_overall_statistics["SNP_calls"]["CL"][row['CELLS']] += local_counter

    print("Made {}/{} VCFs ({} experiment VCFs, {} control VCFs), {} annotated tables".format(
        made_control_vcfs + made_experiment_vcfs, total_vcf_count,
        made_experiment_vcfs,
        made_control_vcfs,
        made_annotated_tables))

    obj_counter = {'TF': 0, 'CL': 0}
    total_fdrs_ref = 0
    total_fdrs_alt = 0
    for what_for in ('CL', ):
        for obj in os.listdir(get_result_dir_path(what_for)):
            obj_counter[what_for] += 1
            obj_table = pd.read_table(get_result_table_path(what_for, os.path.splitext(obj)[0]))
            if obj_table.empty:
                continue
            local_counter = len(obj_table.index)
            local_counter_rs = len(obj_table['ID'].unique())
            fdr_counter_ref = len(obj_table[(obj_table['fdrp_bh_ref'] <= 0.05)].index)
            fdr_counter_alt = len(obj_table[(obj_table['fdrp_bh_alt'] <= 0.05)].index)
            if obj not in dict_overall_statistics["unique_SNPs"][what_for]:
                dict_overall_statistics["unique_SNPs"][what_for][obj] = 0
                dict_overall_statistics["unique_SNPs_rs"][what_for][obj] = 0
            dict_overall_statistics["unique_SNPs"][what_for][obj] += local_counter
            dict_overall_statistics["unique_SNPs_rs"][what_for][obj] += local_counter_rs

            if obj not in dict_overall_statistics["unique_asb"][what_for]:
                dict_overall_statistics["unique_asb"][what_for][obj] = 0
                dict_overall_statistics["unique_asb_rs"][what_for][obj] = 0
            dict_overall_statistics["unique_asb"][what_for][obj] += fdr_counter_ref
            dict_overall_statistics["unique_asb_rs"][what_for][obj] += fdr_counter_alt
            total_fdrs_ref += fdr_counter_ref
            total_fdrs_alt += fdr_counter_alt

        print("Made aggregation for {} {}s".format(obj_counter[what_for], what_for))
        print('In {} aggregation - {},{} ref and alt ASB events respectively'.format(what_for,
                                                                                     total_fdrs_ref, total_fdrs_alt))
    with open(os.path.join(get_release_stats_path(), "overall_statistics.json"), "w") as json_file:
        json.dump(dict_overall_statistics, json_file)

    with open(os.path.join(get_release_stats_path(), 'convert_cell_lines.json'), 'w') as o:
        d_to_write = {}
        master_df = pd.read_table(master_list_path, dtype=dtype_dict)
        for key in master_df['CELLS'].tolist():
            d_to_write[key] = remove_punctuation(key)
        json.dump(d_to_write, o, indent=2)


if __name__ == '__main__':
    main()
