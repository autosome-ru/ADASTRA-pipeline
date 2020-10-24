import os
import numpy as np
from scripts.HELPERS.helpers import pack
from scripts.HELPERS.paths_for_components import results_path


def get_concordance(p_val_ref, p_val_alt, motif_fc, motif_pval_ref, motif_pval_alt):
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


def make_dict_from_data(tf_fasta_path, motif_length):
    # read sarus file and choose best hit
    dict_of_snps = {}
    if os.path.isfile(tf_fasta_path):
        motif_length = int(motif_length)
        with open(tf_fasta_path, 'r') as sarus:
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
                    line = line.strip('\n').split("\t")
                    dict_of_snps[current_snp_id][allele].append({
                        "p": float(line[0]),
                        "orientation": line[2],
                        "pos": int(line[1]) if line[2] == '-' else motif_length - 1 - int(line[1]),
                    })
    return dict_of_snps


def main(tf_name, motif_len):
    sarus_dir = os.path.join(results_path, 'Sarus')
    tf_fasta_path = os.path.join(sarus_dir, tf_name + '')
    dict_of_snps = make_dict_from_data(tf_fasta_path, motif_len)
    adjusted_columns = ['motif_log_pref', 'motif_log_palt', 'motif_fc', 'motif_pos', 'motif_orient', "motif_conc"]
    sarus_table_path = os.path.join(sarus_dir, tf_name + '.tsv')
    tf_table_path = os.path.join(results_path, 'TF_P-values', tf_name + '.tsv')
    with open(tf_table_path, 'r') as table, open(sarus_table_path, 'w') as out:
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
                                       get_concordance(float(line[-2]), float(line[-1]),
                                                       motif_fc,
                                                       dict_of_snps[ID]['ref'][best_idx]['p'],
                                                       dict_of_snps[ID]['alt'][best_idx]['p'])
                                       ]))
