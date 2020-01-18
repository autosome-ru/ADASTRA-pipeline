import sys
import os.path
from statistics import median_grouped
from scipy import stats
import numpy as np
import json
import statsmodels.stats.multitest
import pandas as pd
from collections import OrderedDict

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths_for_components import parameters_path, results_path, tf_dict_path, cl_dict_path
from scripts.HELPERS.helpers import callers_names, unpack, pack, check_if_in_expected_args, expected_args, states


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


def read_weights():
    r = {}
    w = {}
    for fixed_allele in ('ref', 'alt'):
        r[fixed_allele] = {}
        w[fixed_allele] = {}
        for BAD in states:
            precalc_params_path = parameters_path + 'NBweights_{}_BAD={:.1f}.npy'.format(fixed_allele, BAD)
            coefs_array = np.load(precalc_params_path)
            r[fixed_allele][BAD] = coefs_array[:, 0]
            w[fixed_allele][BAD] = coefs_array[:, 1]
    return r, w


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
    interesting_dict = {}
    with open("red.txt") as red:
        for line in red:
            interesting_dict[line.strip()] = "red"
    with open("blue.txt") as blue:
        for line in blue:
            interesting_dict[line.strip()] = "blue"
    common_snps = dict()
    for table in tables:
        if os.path.isfile(table):
            with open(table, 'r') as file:
                for line in file:
                    try:
                        (chr, pos, ID, ref, alt, ref_c, alt_c, repeat, in_callers,
                         BAD, dip_qual, lq, rq, seg_c, sum_cov,
                         p_ref, p_alt) = unpack(line, use_in="Aggregation")
                    except ValueError:
                        continue
                    if p_ref == '.' or ID == '.':
                        continue
                    cov = ref_c + alt_c
                    if ID not in interesting_dict:
                        continue

                    try:
                        common_snps[(chr, pos, ID, ref, alt, repeat)].append(
                            (cov, ref_c, alt_c, in_callers, BAD, dip_qual, lq, rq,
                             seg_c, sum_cov,
                             p_ref, p_alt))
                    except KeyError:
                        common_snps[(chr, pos, ID, ref, alt, repeat)] = [
                            (cov, ref_c, alt_c, in_callers, BAD, dip_qual, lq, rq,
                             seg_c, sum_cov,
                             p_ref, p_alt)]

    print('Writing {}'.format(key_name))

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
        sys.exit(0)
    keys = list(filtered_snps.keys())
    debug_dict = {}
    final_dict = {}
    for tr in range(50, 1401, 50):
        debug_dict[tr] = {}
        final_dict[tr] = {}
        for key in keys:
            chr, pos, ID, ref, alt, repeat = key
            debug_dict[tr][ID] = {"col": interesting_dict[ID]}

            value = filtered_snps[key]
            SNP_counter += 1
            if SNP_counter % 10000 == 0:
                print('done {}'.format(SNP_counter))
            pref_array = []
            palt_array = []
            for v in value:
                cov, ref_c, alt_c, in_callers, BAD, dip_qual, lq, rq, seg_c, sum_cov, p_ref, p_alt = v

                if cov >= tr:
                    continue
                pref_array.append(p_ref)
                palt_array.append(p_alt)

            logitp_ref = logit_combine_p_values(pref_array)
            logitp_alt = logit_combine_p_values(palt_array)
            debug_dict[tr][ID]["p"] = min(logitp_ref, logitp_alt)
        strange_list = set(debug_dict[tr][x]["p"] for x in debug_dict[tr])
        strange_list_len = len(strange_list)
        for strange_tr in strange_list:
            final_dict[tr][strange_tr] = [len([x for x in debug_dict[tr]
                                               if (debug_dict[tr][x]['col'] == 'blue' and
                                                   debug_dict[tr][x]['p'] <= strange_tr)]),
                                          len([x for x in debug_dict[tr] if
                                               (debug_dict[tr][x]['col'] == 'red' and
                                                debug_dict[tr][x]['p'] <= strange_tr)])]
        x = list(strange_list)
        print([final_dict[tr][x_] for x_ in x])
        y1 = [final_dict[tr][x_][0] for x_ in x]
        y2 = [final_dict[tr][x_][1] for x_ in x]
        final_dict[tr] = {}
        final_dict[tr]['x'] = x
        final_dict[tr]['y1'] = y1
        final_dict[tr]['y2'] = y2

    with open(parameters_path + "debug_dict.json", "w") as out:
        json.dump(final_dict, out)
