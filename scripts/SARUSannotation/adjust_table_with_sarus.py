import os
import sys
import numpy as np
from scripts.HELPERS.helpers import pack


def get_color(p_val_ref, p_val_alt, motif_fc, motif_pval_ref, motif_pval_alt):
    log_pv = np.log10(min(p_val_ref, p_val_alt)) * np.sign(p_val_alt - p_val_ref)
    if abs(log_pv) < -np.log10(0.05) or motif_fc == 0:
        return None
    if max(motif_pval_ref, motif_pval_alt) >= -np.log10(0.0005):
        result = "Weak " if abs(motif_fc) < 2 else ""
        if motif_fc * log_pv > 0:
            result += 'Concordant'
        elif motif_fc * log_pv < 0:
            result += 'Discordant'
        return result
    else:
        return "No Hit"


# read sarus file and choose best
dict_of_snps = {}
if os.path.isfile(sys.argv[2]):
    motif_length = int(sys.argv[4])
    with open(sys.argv[2], 'r') as sarus:
        allele = None
        current_snp_id = None
        for line in sarus:
            if line[0] == ">":
                # choose best
                allele = line[-4:-1]
                current_snp_id = line[1:-5]
                assert allele in ("ref", "alt")
                if allele == "ref":
                    dict_of_snps[current_snp_id] = {"ref": [], "alt": []}
            else:
                assert allele in ("ref", "alt")
                line = line.strip().split("\t")
                dict_of_snps[current_snp_id][allele].append({
                    "p": float(line[0]),
                    "orientation": line[2],
                    "pos": int(line[1]) if line[2] == '-' else motif_length - 1 - int(line[1]),
                })

adjusted_columns = ['motif_log_pref', 'motif_log_palt', 'motif_fc', 'motif_pos', 'motif_orient', "motif_conc"]
with open(sys.argv[1], 'r') as table, open(sys.argv[3], 'w') as out:
    for line in table:
        line = line.strip('\n').split('\t')
        if line[0][0] == '#':
            out.write(pack(line + adjusted_columns))
            continue
        if len(dict_of_snps) == 0:
            out.write(pack(line + [""] * len(adjusted_columns)))
            continue
        ID = line[2] + ";" + line[4]

        dict_of_snps[ID]['ref'] = sorted(dict_of_snps[ID]['ref'], key=lambda x: x['pos'])
        dict_of_snps[ID]['ref'] = sorted(dict_of_snps[ID]['ref'], key=lambda x: x['orientation'])

        dict_of_snps[ID]['alt'] = sorted(dict_of_snps[ID]['alt'], key=lambda x: x['pos'])
        dict_of_snps[ID]['alt'] = sorted(dict_of_snps[ID]['alt'], key=lambda x: x['orientation'])

        ref_best = max(enumerate(dict_of_snps[ID]['ref']), key=lambda x: x[1]['p'])
        alt_best = max(enumerate(dict_of_snps[ID]['alt']), key=lambda x: x[1]['p'])

        best_idx, _ = max((ref_best, alt_best), key=lambda x: x[1]['p'])

        if len(dict_of_snps[ID]['ref']) != len(dict_of_snps[ID]['alt']):
            print(ID, dict_of_snps[ID]['ref'], dict_of_snps[ID]['alt'])

        assert len(dict_of_snps[ID]['ref']) == len(dict_of_snps[ID]['alt'])
        assert dict_of_snps[ID]['ref'][best_idx]['pos'] == dict_of_snps[ID]['alt'][best_idx]['pos']
        motif_fc = (dict_of_snps[ID]['alt'][best_idx]['p'] - dict_of_snps[ID]['ref'][best_idx]['p']) / np.log10(2)
        if line[-1] == "":
            out.write(pack(line + [dict_of_snps[ID]['ref'][best_idx]['p'],
                                   dict_of_snps[ID]['alt'][best_idx]['p'],
                                   motif_fc,
                                   dict_of_snps[ID]['ref'][best_idx]['pos'],
                                   dict_of_snps[ID]['ref'][best_idx]['orientation'], ""]))
        else:
            out.write(pack(line + [dict_of_snps[ID]['ref'][best_idx]['p'],
                                   dict_of_snps[ID]['alt'][best_idx]['p'],
                                   motif_fc,
                                   dict_of_snps[ID]['ref'][best_idx]['pos'],
                                   dict_of_snps[ID]['ref'][best_idx]['orientation'],
                                   get_color(float(line[-2]), float(line[-1]),
                                             motif_fc,
                                             dict_of_snps[ID]['ref'][best_idx]['p'],
                                             dict_of_snps[ID]['alt'][best_idx]['p'])
                                   ]))
