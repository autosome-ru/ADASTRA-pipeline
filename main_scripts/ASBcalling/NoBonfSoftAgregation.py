import sys
import os
from statistics import median_grouped
from scipy import stats
import numpy as np


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
    return chr, pos, ID, ref, alt, ref_c, alt_c, Q, GQ, in_macs, in_sissrs, in_cpics, in_gem, callers, ploidy, dip_qual, lq, rq, seg_c, p_ref, p_alt


def pack(values):
    return '\t'.join(map(str, values)) + '\n'


def correctP(p, BonF):
    p = max(p, 0)
    return min(1, p * BonF)


results_path = '/home/abramov/RESULTS/TF_NoCorr_P-values/'




    print('Reading ' + TF + ' k = {}'.format(k))
    common_snps = dict()
    for table in tables:
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
                                                                  seg_c, p_ref, p_alt))
                except KeyError:
                    common_snps[(chr, pos, ID, ref, alt)] = [(cov, ref_c, alt_c, callers, ploidy, dip_qual, lq, rq,
                                                              seg_c, p_ref, p_alt)]

    print('Writing', TF)
    with open(results_path + TF + '_common_table.tsv', 'w') as out:
        out.write(pack(['#chr', 'pos', 'ID', 'ref', 'alt', 'm_callers', 'm_ploidy', 'm_q', 'm_dipq',
                        'm_segc', 'm_datasets', 'maxdepth_ref/alt', 'maxdepth_ploidy', 'maxdepth_m1',
                        'maxdepth_m2', 'mostsig_ref/alt', 'mostsig_ploidy', 'mostsig_m1', 'mostsig_m2',
                        'min_cover', 'max_cover', 'med_cover', 'mean_cover', 'total_cover', 'm1', 'm2',
                        'm_hpref', 'm_hpalt', 'm_fpref', 'm_fpalt', 'm_stpref', 'm_stpalt']))
        print(len(common_snps), 'snps')
        filtered_snps = dict()
        BonF = 1
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

        assert (BonF == 1)
        print(BonF, ' filtered')
        counter = 0
        keys = list(filtered_snps.keys())
        keys = sorted(keys, key=lambda x: x[1])
        keys = sorted(keys, key=lambda x: x[0])
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

            for v in value:
                (cov, ref_c, alt_c, callers, ploidy, dip_qual, lq, rq, seg_c, p_ref, p_alt) = v
                if p_ref <= 0 or p_alt <= 0:
                    print('Wrong P!')
                    continue
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
            m_hpref = correctP(stats.hmean(c_pref), BonF)
            m_hpalt = correctP(stats.hmean(c_palt), BonF)
            m_fpref = correctP(stats.combine_pvalues(c_pref, method='fisher')[1], BonF)
            m_fpalt = correctP(stats.combine_pvalues(c_palt, method='fisher')[1], BonF)
            m_stpref = correctP(stats.combine_pvalues(c_pref, method='stouffer')[1], BonF)
            m_stpalt = correctP(stats.combine_pvalues(c_palt, method='stouffer')[1], BonF)
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
