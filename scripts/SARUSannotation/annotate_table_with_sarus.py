import os
import numpy as np
import pandas as pd
from scripts.HELPERS.paths_for_components import results_path
from scripts.HELPERS.paths import get_tf_sarus_path

cols = ['motif_log_pref',
        'motif_log_palt',
        'motif_fc',
        'motif_pos',
        'motif_orient',
        'motif_conc']


def get_concordance(p_val_ref, p_val_alt, motif_fc, motif_pval_ref, motif_pval_alt):
    if not pd.isna(p_val_ref) and not pd.isna(p_val_alt):
        log_pv = np.log10(min(p_val_ref, p_val_alt)) * np.sign(p_val_alt - p_val_ref)
        if abs(log_pv) >= -np.log10(0.25):
            if max(motif_pval_ref, motif_pval_alt) >= -np.log10(0.0005) and motif_fc != 0:
                result = "Weak " if abs(motif_fc) < 2 else ""
                if motif_fc * log_pv > 0:
                    result += 'Concordant'
                elif motif_fc * log_pv < 0:
                    result += 'Discordant'
                return result
            else:
                return "No Hit"
    return None


def make_dict_from_data(tf_fasta_path, motif_length):
    # read sarus file and choose best hit
    dict_of_snps = {}
    if os.path.isfile(tf_fasta_path) and motif_length is not None:
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
                    line = line.strip('\n').split()
                    dict_of_snps[current_snp_id][allele].append({
                        "p": float(line[0]),
                        "orientation": line[2],
                        "pos": int(line[1]) if line[2] == '-' else motif_length - 1 - int(line[1]),
                    })
    return dict_of_snps


def adjust_with_sarus(df_row, dict_of_snps):
    ID = "{};{}".format(df_row['ID'], df_row['alt'])
    if len(dict_of_snps) == 0:

        result = [None] * 6
    else:
        dict_of_snps[ID]['ref'] = sorted(dict_of_snps[ID]['ref'], key=lambda x: x['pos'])
        dict_of_snps[ID]['ref'] = sorted(dict_of_snps[ID]['ref'], key=lambda x: x['orientation'])
        dict_of_snps[ID]['alt'] = sorted(dict_of_snps[ID]['alt'], key=lambda x: x['pos'])
        dict_of_snps[ID]['alt'] = sorted(dict_of_snps[ID]['alt'], key=lambda x: x['orientation'])
        ref_best = max(enumerate(dict_of_snps[ID]['ref']), key=lambda x: x[1]['p'])
        alt_best = max(enumerate(dict_of_snps[ID]['alt']), key=lambda x: x[1]['p'])
        best_idx, _ = max((ref_best, alt_best), key=lambda x: x[1]['p'])

        if len(dict_of_snps[ID]['ref']) != len(dict_of_snps[ID]['alt']):
            raise AssertionError(ID, dict_of_snps[ID]['ref'], dict_of_snps[ID]['alt'])
        assert dict_of_snps[ID]['ref'][best_idx]['pos'] == dict_of_snps[ID]['alt'][best_idx]['pos']
        result = [dict_of_snps[ID]['ref'][best_idx]['p'],
                  dict_of_snps[ID]['alt'][best_idx]['p'],
                  (dict_of_snps[ID]['alt'][best_idx]['p'] - dict_of_snps[ID]['ref'][best_idx]['p']) / np.log10(2),
                  dict_of_snps[ID]['ref'][best_idx]['pos'],
                  dict_of_snps[ID]['ref'][best_idx]['orientation'],
                  get_concordance(df_row['fdrp_bh_ref'],
                                  df_row['fdrp_bh_alt'],
                                  (dict_of_snps[ID]['alt'][best_idx]['p'] -
                                   dict_of_snps[ID]['ref'][best_idx]['p']) / np.log10(2),
                                  dict_of_snps[ID]['ref'][best_idx]['p'],
                                  dict_of_snps[ID]['alt'][best_idx]['p'])
                  ]
    return pd.Series(dict(zip(cols, result)))


def main(tf_name, motif_length):
    after_sarus_fasta_path = get_tf_sarus_path(tf_name, 'sarus')
    dict_of_snps = make_dict_from_data(after_sarus_fasta_path, motif_length)
    tf_df_path = os.path.join(results_path, 'TF_P-values', tf_name + '.tsv')
    tf_df = pd.read_table(tf_df_path)
    tf_df[cols] = tf_df.apply(lambda x:
                              adjust_with_sarus(x, dict_of_snps), axis=1)
    tf_df.to_csv(tf_df_path, header=True, sep='\t', index=False)
