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
from scripts.HELPERS.paths import results_path, cl_dict_path, tf_dict_path, parameters_path
from scripts.HELPERS.helpers import callers_names, unpack, pack, check_if_in_expected_args, expected_args, states


def logit_combine_p_values(pvalues):
    if 0 in pvalues:
        return 0
    pvalues = np.array([p for p in pvalues if 1 > p > 0])
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
                         ploidy, dip_qual, lq, rq, seg_c, sum_cov,
                         p_ref, p_alt) = unpack(line, use_in="Aggregation")
                    except ValueError:
                        continue
                    if p_ref == '.':
                        continue
                    cov = ref_c + alt_c

                    try:
                        common_snps[(chr, pos, ID, ref, alt, repeat)].append(
                            (cov, ref_c, alt_c, in_callers, ploidy, dip_qual, lq, rq,
                             seg_c, sum_cov,
                             p_ref, p_alt,
                             table_name, another_agr))
                    except KeyError:
                        common_snps[(chr, pos, ID, ref, alt, repeat)] = [
                            (cov, ref_c, alt_c, in_callers, ploidy, dip_qual, lq, rq,
                             seg_c, sum_cov,
                             p_ref, p_alt,
                             table_name, another_agr)]

    print('Writing {}'.format(key_name))

    with open(results_path + what_for + "_P-values/" + key_name + '_common_table.tsv', 'w') as out:
        out.write(pack(['#chr', 'pos', 'ID', 'ref', 'alt', 'repeat_type', 'n_peak_calls', 'n_peak_callers',
                        'mean_sBAD',
                        'mean_deltaL_neighborBAD', 'mean_deltaL_BAD1', 'mean_SNP_per_segment', 'n_aggregated',
                        'refc_maxdepth', 'altc_maxdepth', 'sBAD_maxdepth', 'm1_maxdepth', 'm2_maxdepth',
                        'refc_mostsig', 'altc_mostsig', 'sBAD_mostsig', 'm1_mostsig', 'm2_mostsig',
                        'min_cover', 'max_cover', 'median_cover', 'total_cover',
                        'm1_mean_ref', 'm1_mean_alt',
                        'm2_mean_ref', 'm2_mean_alt',
                        'logitp_ref', 'logitp_alt',
                        'fisherp_ref', 'fisherp_alt']))

        filtered_snps = dict()
        for key in common_snps:
            values = []
            accept = False
            for value in common_snps[key]:
                if value[0] >= 10:
                    values.append(value)
                if value[0] >= 10:
                    accept = True
            if accept or sum(value[0] for value in values) >= 10:
                filtered_snps[key] = values

        counter = 0
        print('{} snps'.format(len(filtered_snps)))

        if len(filtered_snps) == 0:
            sys.exit(0)
        origin_of_snp_dict = OrderedDict()
        keys = list(filtered_snps.keys())
        keys = sorted(keys, key=lambda chr_pos: chr_pos[1])
        keys = sorted(keys, key=lambda chr_pos: chr_pos[0])
        for key in keys:
            chr, pos, ID, ref, alt, repeat = key
            value = filtered_snps[key]
            counter += 1
            if counter % 10000 == 0:
                print('done {}'.format(counter))
            c_uniq_callers = dict(zip(callers_names, [False] * len(callers_names)))
            m_total_callers = 0
            c_ploidy = []
            c_dipq = []
            c_q = []
            c_segc = []
            c_pref = []
            c_palt = []
            c_cover = []
            c_m1 = []
            c_m2 = []
            c_table_names = []
            c_another_agr = []
            c_ref = []
            c_alt = []

            for v in value:
                cov, ref_c, alt_c, in_callers, ploidy, dip_qual, lq, rq, seg_c, sum_cov, p_ref, p_alt, table_name, \
                another_agr = v

                c_table_names.append(table_name)
                c_another_agr.append(another_agr)
                for caller in callers_names:
                    c_uniq_callers[caller] = c_uniq_callers[caller] or in_callers[caller]
                    m_total_callers += in_callers[caller]
                c_ploidy.append(ploidy)
                c_dipq.append(dip_qual)
                if ploidy == 1:
                    c_q.append(rq)
                elif ploidy == 5:
                    c_q.append(lq)
                else:
                    c_q.append(min(lq, rq))
                c_segc.append(seg_c)
                c_pref.append(p_ref)
                c_palt.append(p_alt)
                c_cover.append(cov)

                c_ref.append(ref_c)
                c_alt.append(alt_c)

                x = ref_c / (ref_c + alt_c)

                p = 1 / (ploidy + 1)

                if x <= p or x >= 1 - p:
                    c_m1.append(-1 * np.math.log(min(x, 1 - x) / p, 2) * np.sign(ref_c - alt_c))
                    c_m2.append(np.math.log(max(x, 1 - x) / (1 - p), 2) * np.sign(ref_c - alt_c))

            min_cover = min(c_cover)
            max_cover = max(c_cover)
            med_cover = median_grouped(c_cover)
            total_cover = sum(c_cover)
            m_unique_callers = sum(c_uniq_callers[caller] for caller in callers_names)
            m_ploidy = np.round(np.mean(c_ploidy), 2)
            m_dipq = np.round(np.mean(c_dipq), 1)
            m_q = np.round(np.mean(c_q), 1)
            m_segc = np.round(np.mean(c_segc), 1)
            m_datasets = len(value)

            fisherp_ref = stats.combine_pvalues(c_pref)[1]
            fisherp_alt = stats.combine_pvalues(c_palt)[1]
            m_logpref = logit_combine_p_values(c_pref)
            m_logpalt = logit_combine_p_values(c_palt)

            c_m1_ref = [x for x in c_m1 if x > 0]
            if c_m1_ref:
                m1_ref = np.round(np.mean(c_m1_ref), 3)
            else:
                m1_ref = 0

            c_m1_alt = [x for x in c_m1 if x < 0]
            if c_m1_alt:
                m1_alt = np.round(np.mean(c_m1_alt), 3)
            else:
                m1_alt = 0

            c_m2_ref = [x for x in c_m2 if x > 0]
            if c_m2_ref:
                m2_ref = np.round(np.mean(c_m2_ref), 3)
            else:
                m2_ref = 0

            c_m2_alt = [x for x in c_m2 if x < 0]
            if c_m2_alt:
                m2_alt = np.round(np.mean(c_m2_alt), 3)
            else:
                m2_alt = 0

            m1_dict = dict()
            m2_dict = dict()
            p_dict = dict()
            ref_dict = dict()
            alt_dict = dict()

            for method, sort_key in (('maxdepth', lambda j: -1 * c_cover[j]),
                                     ('mostsig', lambda j: min(c_pref[j], c_palt[j]))):
                try:
                    i_most = min([i for i in range(len(c_cover))
                                  if np.sign(c_ref[i] - c_alt[i]) == np.sign(m_logpalt - m_logpref)],
                                 key=sort_key)
                except ValueError:
                    ref_dict[method] = 'NaN'
                    alt_dict[method] = 'NaN'
                    p_dict[method] = 'NaN'
                    m1_dict[method] = 'NaN'
                    m2_dict[method] = 'NaN'
                    continue
                ref_dict[method] = str(c_ref[i_most])
                alt_dict[method] = str(c_alt[i_most])
                p_dict[method] = c_ploidy[i_most]
                x = c_ref[i_most] / (c_ref[i_most] + c_alt[i_most])
                p = 1 / (c_ploidy[i_most] + 1)
                if x <= p or x >= 1 - p:
                    m1_dict[method] = -1 * np.math.log(min(x, 1 - x) / p, 2) * np.sign(c_ref[i_most] - c_alt[i_most])
                    m2_dict[method] = np.math.log(max(x, 1 - x) / (1 - p), 2) * np.sign(c_ref[i_most] - c_alt[i_most])
                else:
                    m1_dict[method] = 0
                    m2_dict[method] = 0

            out.write(pack(
                [chr, pos, ID, ref, alt, repeat, m_total_callers, m_unique_callers,
                 m_ploidy, m_q, m_dipq, m_segc, m_datasets,
                 ref_dict['maxdepth'], alt_dict['maxdepth'], p_dict['maxdepth'],
                 m1_dict['maxdepth'], m2_dict['maxdepth'],
                 ref_dict['mostsig'], alt_dict['mostsig'], p_dict['mostsig'],
                 m1_dict['mostsig'], m2_dict['mostsig'],
                 min_cover, max_cover, med_cover, total_cover,
                 m1_ref, abs(m1_alt), m2_ref, abs(m2_alt),
                 m_logpref, m_logpalt,
                 fisherp_ref, fisherp_alt]))
            origin_of_snp_dict["\t".join(map(str, key))] = {'aligns': c_table_names,
                                                            expected_args[what_for]: c_another_agr,
                                                            'ref_counts': c_ref, 'alt_counts': c_alt,
                                                            'ref_pvalues': c_pref, 'alt_pvalues': c_palt}

    print("Counting FDR")
    table = pd.read_table(results_path + what_for + "_P-values/" + key_name + '_common_table.tsv')
    if table.empty:
        sys.exit(0)

    bool_ar_ref, p_val_ref, _, _ = statsmodels.stats.multitest.multipletests(table["logitp_ref"],
                                                                             alpha=0.05, method='fdr_bh')
    bool_ar_alt, p_val_alt, _, _ = statsmodels.stats.multitest.multipletests(table["logitp_alt"],
                                                                             alpha=0.05, method='fdr_bh')
    table["fdrp_ref"] = pd.Series(p_val_ref)
    table["fdrp_alt"] = pd.Series(p_val_alt)

    bool_ar_ref_fisher, p_val_ref_fisher, _, _ = statsmodels.stats.multitest.multipletests(table["fisherp_ref"],
                                                                                           alpha=0.05, method='fdr_bh')
    bool_ar_alt_fisher, p_val_alt_fisher, _, _ = statsmodels.stats.multitest.multipletests(table["fisherp_alt"],
                                                                                           alpha=0.05, method='fdr_bh')
    table["fdrp_ref_fisher"] = pd.Series(p_val_ref_fisher)
    table["fdrp_alt_fisher"] = pd.Series(p_val_alt_fisher)

    with open(results_path + what_for + "_P-values/" + key_name + '_common_table.tsv', "w") as w:
        table.to_csv(w, sep="\t", index=False)

    bool_ar = bool_ar_ref + bool_ar_alt

    annotate_snp_with_tables(origin_of_snp_dict, p_val_ref, p_val_alt, bool_ar)

    with open(results_path + what_for + '_DICTS/' + key_name + '_DICT.json', 'w') as out:
        json.dump(origin_of_snp_dict, out)
