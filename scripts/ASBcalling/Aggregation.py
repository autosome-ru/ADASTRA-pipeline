import sys
import os.path
from statistics import median_grouped
from scipy import stats
import numpy as np
import json
import statsmodels.stats.multitest
import pandas as pd
from collections import OrderedDict
from scripts.HELPERS.paths_for_components import results_path, tf_dict_path, cl_dict_path
from scripts.HELPERS.helpers import callers_names, unpack, pack, check_if_in_expected_args, \
    expected_args


def logit_combine_p_values(pvalues):
    pvalues = np.array([pvalue for pvalue in pvalues if 1 > pvalue > 0])
    if len(pvalues) == 0:
        return 1
    elif len(pvalues) == 1:
        return pvalues[0]

    statistic = -np.sum(np.log(pvalues)) + np.sum(np.log1p(-pvalues))
    k = len(pvalues)
    nu = np.int_(5 * k + 4)
    approx_factor = np.sqrt(np.int_(3) * nu / (np.int_(k) * np.square(np.pi) * (nu - np.int_(2))))
    pval = stats.distributions.t.sf(statistic * approx_factor, nu)
    return pval


def annotate_snp_with_tables(dictionary, ps_ref, ps_alt, bool_ar):  # return part of the dictionary with fdr from table
    keys = list(dictionary.keys())
    for index in range(len(ps_ref)):
        key = keys[index]
        if bool_ar[index]:
            dictionary[key]['logitp_ref'] = ps_ref[index]
            dictionary[key]['logitp_alt'] = ps_alt[index]
        else:
            del dictionary[key]


def get_name(path):  # path format */ALIGNS000000_table_p.txt
    return path.split("/")[-1].split("_")[0]


def invert(dictionary):
    inverted_dictionary = {}
    for key in dictionary:
        for value in dictionary[key]:
            inverted_dictionary[value] = key
    return inverted_dictionary


def get_another_agr(path, what_for):
    if what_for == "TF":
        return invert(cell_lines_dict).get(path, "None")
    if what_for == "CL":
        return invert(tf_dict).get(path, "None")


def get_noise(k, n, weight):
    if n == 10:
        return 1
    if k <= 4 or k >= n - 4:
        return 0
    return weight * max(k - n / 2, 0) / (1 / 2 * (n - 5 - n // 2) * (n // 2 - 4))


if __name__ == '__main__':

    what_for = sys.argv[1]  # "TF" or "CL" arguments are expected
    check_if_in_expected_args(what_for)
    key_name = sys.argv[2]

    if not os.path.isdir(results_path + what_for + '_DICTS/'):
        os.mkdir(results_path + what_for + '_DICTS/')
    if not os.path.isdir(results_path + what_for + "_P-values/"):
        os.mkdir(results_path + what_for + "_P-values/")

    table_path = results_path + what_for + '_P-values/{}.tsv'.format(key_name)
    with open(cl_dict_path, "r") as read_file:
        cell_lines_dict = json.loads(read_file.readline())
    with open(tf_dict_path, "r") as read_file:
        tf_dict = json.loads(read_file.readline())
    tables = []
    if what_for == "CL":
        tables = cell_lines_dict[key_name]
    if what_for == "TF":
        tables = tf_dict[key_name]
    print('Reading datasets for {} {}'.format(what_for, key_name))
    common_snps = dict()
    for table in tables:
        if os.path.isfile(table):
            table_name = get_name(table)
            another_agr = get_another_agr(table,
                                          what_for)  # returns name of cell-line for aggregation on TF and vice versa
            with open(table, 'r') as file:
                for line in file:
                    try:
                        (chr, pos, ID, ref, alt, ref_c, alt_c, repeat, in_callers,
                         BAD, Quals, seg_c, sum_cov, p_ref, p_alt, es_ref, es_alt) = unpack(line, use_in="Aggregation")
                    except ValueError:
                        if line.startswith('#'):
                            continue
                        else:
                            raise
                    if p_ref == '.' or ID == '.':
                        continue
                    cov = ref_c + alt_c

                    try:
                        common_snps[(chr, pos, ID, ref, alt, repeat)].append(
                            (cov, ref_c, alt_c, in_callers, BAD, Quals,
                             seg_c, sum_cov,
                             p_ref, p_alt, es_ref, es_alt,
                             table_name, another_agr))
                    except KeyError:
                        common_snps[(chr, pos, ID, ref, alt, repeat)] = [
                            (cov, ref_c, alt_c, in_callers, BAD, Quals,
                             seg_c, sum_cov,
                             p_ref, p_alt, es_ref, es_alt,
                             table_name, another_agr)]
        else:
            print("There is no {}".format(table))
    print('Writing {}'.format(key_name))

    with open(table_path, 'w') as out:
        out.write(pack(['#chr', 'pos', 'ID', 'ref', 'alt', 'repeat_type', 'n_peak_calls', 'n_peak_callers',
                        'mean_BAD', 'mean_SNP_per_segment', 'n_aggregated',
                        'refc_mostsig_ref', 'altc_mostsig_ref', 'BAD_mostsig_ref', 'es_mostsig_ref',
                        'refc_mostsig_alt', 'altc_mostsig_alt', 'BAD_mostsig_alt', 'es_mostsig_alt',
                        'min_cover', 'max_cover', 'median_cover', 'total_cover',
                        'es_mean_ref', 'es_mean_alt',
                        'logitp_ref', 'logitp_alt']))
        filtered_snps = dict()
        for key in common_snps:
            values = []
            accept = False
            for value in common_snps[key]:
                # filtering part is now redundant
                if value[0] >= 10:
                    values.append(value)
                if value[0] >= 10:
                    accept = True
            if accept or sum(value[0] for value in values) >= 10:
                filtered_snps[key] = values

        SNP_counter = 0
        print('{} snps'.format(len(filtered_snps)))

        if len(filtered_snps) == 0:
            os.remove(table_path)
            sys.exit(0)
        origin_of_snp_dict = OrderedDict()
        keys = list(filtered_snps.keys())
        keys = sorted(keys, key=lambda chr_pos: chr_pos[1])
        keys = sorted(keys, key=lambda chr_pos: chr_pos[0])
        for key in keys:
            chr, pos, ID, ref, alt, repeat = key
            value = filtered_snps[key]
            SNP_counter += 1
            if SNP_counter % 10000 == 0:
                print('done {}'.format(SNP_counter))
            uniq_callers_counter = dict(zip(callers_names, [False] * len(callers_names)))
            total_callers_counter = 0
            BAD_array = []
            SNPs_per_segment_array = []
            pref_array = []
            palt_array = []
            cover_array = []
            ref_effect_size_array = []
            alt_effect_size_array = []
            table_names_array = []
            another_agr_name = []
            ref_counts_array = []
            alt_counts_array = []

            for v in value:
                cov, ref_c, alt_c, in_callers, BAD, Quals, seg_c, sum_cov, p_ref, p_alt, es_ref, es_alt, table_name, \
                another_agr = v

                table_names_array.append(table_name)
                another_agr_name.append(another_agr)
                for caller in callers_names:
                    uniq_callers_counter[caller] = uniq_callers_counter[caller] or in_callers[caller]
                    total_callers_counter += in_callers[caller]
                BAD_array.append(BAD)
                SNPs_per_segment_array.append(seg_c)
                pref_array.append(p_ref)
                palt_array.append(p_alt)
                if es_ref is not None:
                    ref_effect_size_array.append(es_ref/np.log(2))
                if es_alt is not None:
                    alt_effect_size_array.append(es_alt/np.log(2))
                cover_array.append(cov)

                ref_counts_array.append(ref_c)
                alt_counts_array.append(alt_c)
                p = 1 / (BAD + 1)

            min_cover = min(cover_array)
            max_cover = max(cover_array)
            med_cover = median_grouped(cover_array)
            total_cover = sum(cover_array)
            unique_callers = sum(uniq_callers_counter[caller] for caller in callers_names)
            mean_BAD = np.round(np.mean(BAD_array), 2)
            mean_SNPs_per_segment = np.round(np.mean(SNPs_per_segment_array), 1)
            n_aggregated = len(value)

            logitp_ref = logit_combine_p_values(pref_array)
            logitp_palt = logit_combine_p_values(palt_array)

            if ref_effect_size_array:
                weights = [-1 * np.log10(x) for x in pref_array if x != 1]
                es_mean_ref = np.round(np.average(ref_effect_size_array, weights=weights), 3)
                es_mostsig_ref = ref_effect_size_array[np.argmax(weights)]
                idx = np.argmax([-x for x in pref_array])
                ref_c_mostsig_ref = ref_counts_array[idx]
                alt_c_mostsig_ref = alt_counts_array[idx]
                BAD_mostsig_ref = BAD_array[idx]
            else:
                es_mean_ref = 'NaN'
                es_mostsig_ref = 'NaN'
                ref_c_mostsig_ref = 'NaN'
                alt_c_mostsig_ref = 'NaN'
                BAD_mostsig_ref = 'NaN'

            if alt_effect_size_array:
                weights = [-1 * np.log10(x) for x in palt_array if x != 1]
                es_mean_alt = np.round(np.average(alt_effect_size_array, weights=weights), 3)
                es_mostsig_alt = alt_effect_size_array[np.argmax(weights)]
                idx = np.argmax([-x for x in palt_array])
                ref_c_mostsig_alt = ref_counts_array[idx]
                alt_c_mostsig_alt = alt_counts_array[idx]
                BAD_mostsig_alt = BAD_array[idx]
            else:
                es_mean_alt = 'NaN'
                es_mostsig_alt = 'NaN'
                ref_c_mostsig_alt = 'NaN'
                alt_c_mostsig_alt = 'NaN'
                BAD_mostsig_alt = 'NaN'

            out.write(pack(
                [chr, pos, ID, ref, alt, repeat, total_callers_counter, unique_callers,
                 mean_BAD, mean_SNPs_per_segment, n_aggregated,
                 ref_c_mostsig_ref, alt_c_mostsig_ref, BAD_mostsig_ref, es_mostsig_ref,
                 ref_c_mostsig_alt, alt_c_mostsig_alt, BAD_mostsig_alt, es_mostsig_alt,
                 min_cover, max_cover, med_cover, total_cover,
                 es_mean_ref, es_mean_alt,
                 logitp_ref, logitp_palt]))
            origin_of_snp_dict["\t".join(map(str, key))] = {'aligns': table_names_array,
                                                            expected_args[what_for]: another_agr_name,
                                                            'ref_counts': ref_counts_array,
                                                            'alt_counts': alt_counts_array,
                                                            'ref_ef': ref_effect_size_array,
                                                            'alt_ef': alt_effect_size_array,
                                                            'BAD': BAD_array,
                                                            'ref_pvalues': pref_array,
                                                            'alt_pvalues': palt_array,
                                                            }

    print("Counting FDR")

    table = pd.read_table(table_path)
    if table.empty:
        os.remove(table_path)
        sys.exit(0)

    mc_filter_array = np.array(table['max_cover'] >= 20)
    if sum(mc_filter_array) != 0:
        bool_ar_ref, p_val_ref, _, _ = statsmodels.stats.multitest.multipletests(table[mc_filter_array]["logitp_ref"],
                                                                                 alpha=0.05, method='fdr_bh')
        bool_ar_alt, p_val_alt, _, _ = statsmodels.stats.multitest.multipletests(table[mc_filter_array]["logitp_alt"],
                                                                                 alpha=0.05, method='fdr_bh')
    else:
        p_val_ref = []
        p_val_alt = []
        bool_ar_ref = []
        bool_ar_alt = []

    fdr_by_ref = np.array(['NaN'] * len(table.index), dtype=np.float128)
    fdr_by_ref[mc_filter_array] = p_val_ref
    table["fdrp_bh_ref"] = fdr_by_ref

    fdr_by_alt = np.array(['NaN'] * len(table.index), dtype=np.float128)
    fdr_by_alt[mc_filter_array] = p_val_alt
    table["fdrp_bh_alt"] = fdr_by_alt

    table.to_csv(table_path, sep="\t", index=False)

    bool_ar = np.array([False] * len(table.index), dtype=np.bool)
    bool_ar[mc_filter_array] = bool_ar_alt + bool_ar_ref

    with open(results_path + what_for + '_DICTS/{}_DICT.json'.format(key_name), 'w') as out:
        json.dump(origin_of_snp_dict, out)
