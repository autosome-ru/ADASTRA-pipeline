import pandas as pd
import subprocess
from scripts.HELPERS.paths_for_components import FA
from scripts.HELPERS.paths import get_tf_sarus_path, get_result_table_path


def get_name(row):
    return [f" {row['ID']}@{row['alt']}@{allele}" for allele in ('ref', 'alt')]


def tf_to_bed(df, motif_length):
    df['start'] = df['pos'] - motif_length
    df['end'] = df['pos'] + motif_length - 1
    df['name'] = df.apply(get_name, axis=1)
    df = df.explode('name')
    df['nuc'] = df.apply(lambda x: x['alt'] if x['name'].endswith('alt') else x['ref'], axis=1)
    return df[['#chr', 'start', 'end', 'name', 'nuc']]


def main(tf_name, motif_length, opened_df=None):
    if motif_length is None:
        exit(0)
    tf_path = get_result_table_path('TF', tf_name)
    if opened_df is not None:
        tf_df = opened_df
    else:
        tf_df = pd.read_table(tf_path)
    if tf_df.empty:
        exit(0)
    bed_df = tf_to_bed(tf_df, motif_length)
    id_nuc_cor = pd.Series(bed_df['nuc'].values, index=bed_df['name'].values).to_dict()
    bed_path = get_tf_sarus_path(tf_name, 'tsv')
    bed_df[['#chr', 'start', 'end', 'name']].to_csv(bed_path, sep='\t', index=None)
    out_path = get_tf_sarus_path(tf_name, 'fasta')
    fasta_buf = subprocess.check_output(['bedtools', 'getfasta', '-fi', FA, '-bed', bed_path, '-name']).decode('utf-8')
    key = None
    with open(out_path, 'w') as out:
        for line in fasta_buf.split('\n'):
            if line.startswith('>'):
                key = line.strip()[1:]
            elif key:
                assert len(line) == motif_length * 2 - 1
                middle_index = motif_length - 1
                if key.endswith('ref'):
                    assert id_nuc_cor[key] == line[middle_index].upper()
                line = line[:middle_index] + id_nuc_cor[key] + line[middle_index + 1:]
                key = None
            out.write(line + '\n')
