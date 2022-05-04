import pandas as pd
import subprocess
from scripts.HELPERS.paths_for_components import FA
from scripts.HELPERS.paths import get_tf_sarus_path, get_result_table_path


def get_name(row):
    return f" {row['ID']}@{row['alt']}"


def tf_to_bed(tf_df):
    tf_df['start'] = tf_df['pos'] - 1
    tf_df['end'] = tf_df['pos']
    tf_df['name'] = tf_df.apply(get_name, axis=1)
    return tf_df[['#chr', 'start', 'end', 'name']]


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
    bed_df = tf_to_bed(tf_df)
    bed_path = get_tf_sarus_path(tf_name, 'tsv')
    bed_df.to_csv(bed_path, sep='\t', index=None)
    out_path = get_tf_sarus_path(tf_name, 'fasta')
    fasta_buf = subprocess.check_output(['bedtools', 'getfasta', '-fi', FA, '-bed', bed_path, '-name'])
    with open(out_path, 'wb') as out:
        out.write(fasta_buf)
