import sys
import os.path
from statistics import median_grouped
from scipy import stats
import numpy as np
import json
import statsmodels.stats.multitest
import pandas as pd


def unpack(line):
    line = line.split()
    chr = line[0]
    pos = int(line[1])
    ID = line[2]
    ref = line[3]
    alt = line[4]
    Q = float(line[7])
    ref_c, alt_c, GQ, in_macs, in_sissrs, in_cpics, in_gem = map(int, line[5:7] + line[8:13])
    callers = in_macs + in_sissrs + in_cpics + in_gem
    ploidy, dip_qual, lq, rq, seg_c = map(int, line[14:19])
    if line[19] == '.':
        p_ref = '.'
        p_alt = '.'
    else:
        p_ref, p_alt = map(float, line[19:21])

    return chr, pos, ID, ref, alt, ref_c, alt_c, Q, GQ, in_macs, in_sissrs, in_cpics, in_gem, callers, ploidy,\
           dip_qual, lq, rq, seg_c, p_ref, p_alt


def pack(values):
    return '\t'.join(map(str, values)) + '\n'


def annotate_snp_with_tables(dictionary, pandas_table):  # return part of the dictionary with fdr from table
    pass


def get_name(path):  # path format */ALIGNS000000_table_p.txt
    return path.split("/")[-1].split["_"][0]


results_path = '/home/abramov/TF_P-values/'
parameters_path = '/home/abramov/PARAMETERS/'
dicts_path = '/home/abramov/DATA/DICTS/'

what_for = sys.argv[1]
key = sys.argv[2]
with open(parameters_path + what_for + "_DICT.json", "r") as read_file:
    d = json.loads(read_file.readline())  # read CL or TF json
tables = d.get(key, None)
print('Reading datasets for {} '.format(what_for) + key)
common_snps = dict()
for table in tables:
    if os.path.isfile(table):
        table_name = get_name(table)
        with open(table, 'r') as file:
            for line in file:
                if line[0] == '#':
                    continue
                (chr, pos, ID, ref, alt, ref_c, alt_c, Q, GQ, in_macs, in_sissrs, in_cpics, in_gem, callers,
                 ploidy, dip_qual, lq, rq, seg_c, p_ref, p_alt) = unpack(line)
                if p_ref == '.' or ploidy == 0 or ID == '.':
                    continue
                cov = ref_c + alt_c
                try:
                    common_snps[(chr, pos, ID, ref, alt)].append((cov, ref_c, alt_c, callers, ploidy, dip_qual, lq, rq,
                                                                  seg_c, p_ref, p_alt, table_name))
                except KeyError:
                    common_snps[(chr, pos, ID, ref, alt)] = [(cov, ref_c, alt_c, callers, ploidy, dip_qual, lq, rq,
                                                              seg_c, p_ref, p_alt, table_name)]

print('Writing ', key)

with open(results_path + key + '_common_table.tsv', 'w') as out:
    out.write(pack(['#chr', 'pos', 'ID', 'ref', 'alt', 'm_callers', 'm_ploidy', 'm_q', 'm_dipq',
                    'm_segc', 'm_datasets', 'maxdepth_ref/alt', 'maxdepth_ploidy', 'maxdepth_m1',
                    'maxdepth_m2', 'mostsig_ref/alt', 'mostsig_ploidy', 'mostsig_m1', 'mostsig_m2',
                    'min_cover', 'max_cover', 'med_cover', 'mean_cover', 'total_cover', 'm1', 'm2',
                    'm_hpref', 'm_hpalt', 'm_fpref', 'm_fpalt', 'm_stpref', 'm_stpalt']))
    print(len(common_snps), 'snps')
    filtered_snps = dict()
    for key in common_snps:
        is_greater = False
        SNPs_values = []
        for value in common_snps[key]:
            if value[1] < 3 or value[2] < 3:
                continue
            SNPs_values.append(value)
            if value[0] >= 10:
                is_greater = True
        if is_greater:
            filtered_snps[key] = SNPs_values

    counter = 0

    origin_of_snp_dict = {}
    keys = list(filtered_snps.keys())
    keys = sorted(keys, key=lambda chr_pos: chr_pos[1])
    keys = sorted(keys, key=lambda chr_pos: chr_pos[0])
    for key in keys:
        chr, pos, ID, ref, alt = key
        value = filtered_snps[key]
        counter += 1
        if counter % 10000 == 0:
            print(counter, 'done')
        c_callers = []
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

        for v in value:
            (cov, ref_c, alt_c, callers, ploidy, dip_qual, lq, rq, seg_c, p_ref, p_alt, table_name) = v
            if p_ref <= 0 or p_alt <= 0:
                print('Wrong P!')
                continue
            c_table_names.append(table_name)
            c_callers.append(callers)
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

            x = ref_c / (ref_c + alt_c)

            p = 1 / (ploidy + 1)

            if x <= p or x >= 1 - p:
                c_m1.append(-1 * np.math.log(min(x, 1 - x) / p, 2) * np.sign(ref_c - alt_c))
                c_m2.append(np.math.log(max(x, 1 - x) / (1 - p), 2) * np.sign(ref_c - alt_c))

        v = min(value, key=lambda x: min(x[-2], x[-1]))
        mostsig_refalt = str(v[1]) + '/' + str(v[2])
        mostsig_p = str(v[4])
        x = v[1] / (v[1] + v[2])
        p = 1 / (v[4] + 1)
        if x <= p or x >= 1 - p:
            mostsig_m1 = -1 * np.math.log(min(x, 1 - x) / p, 2) * np.sign(v[1] - v[2])
            mostsig_m2 = np.math.log(max(x, 1 - x) / (1 - p), 2) * np.sign(v[1] - v[2])
        else:
            mostsig_m1 = 0
            mostsig_m2 = 0

        v = min(value, key=lambda x: x[0])
        maxdepth_refalt = str(v[1]) + '/' + str(v[2])
        maxdepth_p = str(v[4])
        x = v[1] / (v[1] + v[2])
        p = 1 / (v[4] + 1)
        if x <= p or x >= 1 - p:
            maxdepth_m1 = -1 * np.math.log(min(x, 1 - x) / p, 2) * np.sign(v[1] - v[2])
            maxdepth_m2 = np.math.log(max(x, 1 - x) / (1 - p), 2) * np.sign(v[1] - v[2])
        else:
            maxdepth_m1 = 0
            maxdepth_m2 = 0

        min_cover = min(c_cover)
        max_cover = max(c_cover)
        med_cover = median_grouped(c_cover)
        mean_cover = np.round(np.mean(c_cover), 1)
        m_callers = np.round(np.mean(c_callers), 2)
        m_ploidy = np.round(np.mean(c_ploidy), 2)
        m_dipq = np.round(np.mean(c_dipq), 1)
        m_q = np.round(np.mean(c_q), 1)
        m_segc = np.round(np.mean(c_segc), 1)
        m_datasets = len(value)
        m_hpref = stats.hmean(c_pref)
        m_hpalt = stats.hmean(c_palt)
        m_fpref = stats.combine_pvalues(c_pref, method='fisher')[1]
        m_fpalt = stats.combine_pvalues(c_palt, method='fisher')[1]
        m_stpref = stats.combine_pvalues(c_pref, method='stouffer')[1]
        m_stpalt = stats.combine_pvalues(c_palt, method='stouffer')[1]

        if c_m1 and c_m2:
            m1 = np.round(np.mean(c_m1), 3)
            m2 = np.round(np.mean(c_m2), 3)
        else:
            m1 = 0
            m2 = 0

        out.write(pack(
            [chr, pos, ID, ref, alt, m_callers, m_ploidy, m_q, m_dipq, m_segc, m_datasets, maxdepth_refalt,
             maxdepth_p, maxdepth_m1, maxdepth_m2, mostsig_refalt, mostsig_p, mostsig_m1, mostsig_m2, min_cover,
             max_cover, med_cover, mean_cover, mean_cover * m_datasets, m1, m2, m_hpref, m_hpalt, m_fpref, m_fpalt,
             m_stpref, m_stpalt]))
        origin_of_snp_dict[",".join(map(str, key))] = c_table_names
out.close()
print("Counting FDR")
with open(results_path + key + '_common_table.tsv', 'r') as f:
    table = pd.read_table(f)
    f.close()
    table["m_fdr"] = table[["m_fpref", "m_fpalt"]].min(axis=1)
    bool_ar, p_val = statsmodels.stats.multitest.multipletests(table["m_fdr"],
                                                               alpha=0.05, method='fdr_bh')
    table["m_fdr"] = pd.Series(p_val)
    print(bool_ar)
    with open(results_path + key + '_common_table.tsv', "w") as w:
        table.to_csv(w, sep="\t", index=False)
    table = table.loc(bool_ar)
    datasets_for_SNPs = annotate_snp_with_tables(origin_of_snp_dict, table)
    with open(dicts_path + key + '_DICT.json', 'w') as out:
        json.dump(datasets_for_SNPs, out)
