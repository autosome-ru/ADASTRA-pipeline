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
                    if p_ref == '.' or ID == '.':
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
                        'mean_BAD',
                        'mean_deltaL_neighborBAD', 'mean_deltaL_BAD1', 'mean_SNP_per_segment', 'n_aggregated',
                        'refc_mostsig_ref', 'altc_mostsig_ref', 'BAD_mostsig_ref', 'm_mostsig_ref',
                        'refc_mostsig_alt', 'altc_mostsig_alt', 'BAD_mostsig_alt', 'm_mostsig_alt',
                        'min_cover', 'max_cover', 'median_cover', 'total_cover',
                        'm_mean_ref', 'm_mean_alt',
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
        r, w = read_weights()
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
            c_m_ref = []
            c_m_alt = []
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

                p = 1 / (ploidy + 1)

                if p_ref != 1:
                    if alt_c > 500:
                        r_ref = alt_c
                        w_ref = 0.5
                    else:
                        r_ref = r['alt'][ploidy][alt_c]
                        w_ref = w['alt'][ploidy][alt_c]
                    dist1 = stats.nbinom(r_ref, p)
                    dist2 = stats.nbinom(r_ref, 1 - p)
                    cdf1_ref = dist1.cdf
                    cdf2_ref = dist2.cdf
                    cdf_ref = lambda x: w_ref * cdf1_ref(x) + (1 - w_ref) * cdf2_ref(x)
                    pmf1_ref = dist1.pmf
                    pmf2_ref = dist2.pmf
                    pmf_ref = lambda x: w_ref * pmf1_ref(x) + (1 - w_ref) * pmf2_ref(x)

                    E_ref = (r_ref * (ploidy * w_ref + (1 - w_ref) / ploidy) - sum(i * pmf_ref(i) for i in range(5))) / (
                            1 - cdf_ref(4))

                    c_m_ref.append(np.math.log(ref_c / E_ref))

                if p_alt != 1:
                    if ref_c > 500:
                        r_alt = ref_c
                        w_alt = 0.5
                    else:
                        r_alt = r['ref'][ploidy][ref_c]
                        w_alt = w['ref'][ploidy][ref_c]
                    dist1 = stats.nbinom(r_alt, p)
                    dist2 = stats.nbinom(r_alt, 1 - p)
                    cdf1_alt = dist1.cdf
                    cdf2_alt = dist2.cdf
                    cdf_alt = lambda x: w_alt * cdf1_alt(x) + (1 - w_alt) * cdf2_alt(x)
                    pmf1_alt = dist1.pmf
                    pmf2_alt = dist2.pmf
                    pmf_alt = lambda x: w_alt * pmf1_alt(x) + (1 - w_alt) * pmf2_alt(x)

                    E_alt = (r_alt * (ploidy * w_alt + (1 - w_alt) / ploidy) - sum(i * pmf_alt(i) for i in range(5))) / (
                                1 - cdf_alt(4))

                    c_m_alt.append(np.math.log(alt_c / E_alt))

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

            if c_m_ref:
                weights = [-1 * np.log10(x) for x in c_pref if x != 1]
                m_ref = np.round(np.average(c_m_ref, weights=weights), 3)
                m_mostsig_ref = c_m_ref[np.argmax(weights)]
                idx = np.argmax([-x for x in c_pref])
                ref_c_mostsig_ref = c_ref[idx]
                alt_c_mostsig_ref = c_alt[idx]
                BAD_mostsig_ref = c_ploidy[idx]
            else:
                m_ref = 'NaN'
                m_mostsig_ref = 'NaN'
                ref_c_mostsig_ref = 'NaN'
                alt_c_mostsig_ref = 'NaN'
                BAD_mostsig_ref = 'NaN'

            if c_m_alt:
                weights = [-1 * np.log10(x) for x in c_palt if x != 1]
                m_alt = np.round(np.average(c_m_alt, weights=weights), 3)
                m_mostsig_alt = c_m_alt[np.argmax(weights)]
                idx = np.argmax([-x for x in c_palt])
                ref_c_mostsig_alt = c_ref[idx]
                alt_c_mostsig_alt = c_alt[idx]
                BAD_mostsig_alt = c_ploidy[idx]
            else:
                m_alt = 'NaN'
                m_mostsig_alt = 'NaN'
                ref_c_mostsig_alt = 'NaN'
                alt_c_mostsig_alt = 'NaN'
                BAD_mostsig_alt = 'NaN'

            out.write(pack(
                [chr, pos, ID, ref, alt, repeat, m_total_callers, m_unique_callers,
                 m_ploidy, m_q, m_dipq, m_segc, m_datasets,
                 ref_c_mostsig_ref, alt_c_mostsig_ref, BAD_mostsig_ref, m_mostsig_ref,
                 ref_c_mostsig_alt, alt_c_mostsig_alt, BAD_mostsig_alt, m_mostsig_alt,
                 min_cover, max_cover, med_cover, total_cover,
                 m_ref, m_alt,
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

    mc_filter_array = np.array(table['max_cover'] >= 30)
    if sum(mc_filter_array) != 0:
        bool_ar_ref, p_val_ref, _, _ = statsmodels.stats.multitest.multipletests(table[mc_filter_array]["logitp_ref"],
                                                                                 alpha=0.05, method='fdr_by')
        bool_ar_alt, p_val_alt, _, _ = statsmodels.stats.multitest.multipletests(table[mc_filter_array]["logitp_alt"],
                                                                                 alpha=0.05, method='fdr_by')
    else:
        p_val_ref = []
        p_val_alt = []
        bool_ar_ref = []
        bool_ar_alt = []

    fdr_by_ref = np.array(['NaN'] * len(table.index), dtype=np.float128)
    fdr_by_ref[mc_filter_array] = p_val_ref
    table["fdrp_by_ref"] = fdr_by_ref

    fdr_by_alt = np.array(['NaN'] * len(table.index), dtype=np.float128)
    fdr_by_alt[mc_filter_array] = p_val_alt
    table["fdrp_by_alt"] = fdr_by_alt

    with open(results_path + what_for + "_P-values/" + key_name + '_common_table.tsv', "w") as w:
        table.to_csv(w, sep="\t", index=False)

    bool_ar = np.array([False] * len(table.index), dtype=np.bool)
    bool_ar[mc_filter_array] = bool_ar_alt + bool_ar_ref

    annotate_snp_with_tables(origin_of_snp_dict, p_val_ref, p_val_alt, bool_ar)

    with open(results_path + what_for + '_DICTS/' + key_name + '_DICT.json', 'w') as out:
        json.dump(origin_of_snp_dict, out)
