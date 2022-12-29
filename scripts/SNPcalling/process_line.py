import os
import subprocess

from scripts.HELPERS.paths_for_components import alignments_path


dnase_bams_path = '/mnt/NAS/home/abramov/raw_alignments.GTRD/faire'


def make_bams_list(gtrd_id):
    return [file for file in os.listdir(dnase_bams_path)
            if os.path.splitext(file)[0] == gtrd_id and os.path.splitext(file)[1] == 'bam']


def process_bam(bam, out_path):
    cmd = ['bash', '/home/abramov/faire/ADASTRA-pipeline/scripts/SNPcalling/SNPcalling.sh',
           '-Exp', os.path.join(dnase_bams_path, bam),
           '-Out', out_path]
    print(cmd)
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
    if error:
        print(error.decode('utf-8'))


def main(master_line):
    exp_name, align_name, read_groups, _ = master_line.split('\t')
    print('Processing', align_name)
    out_path = os.path.join(alignments_path)
    downloaded_bams = make_bams_list(align_name)
    if downloaded_bams is None or len(downloaded_bams) < 2:
        return

    for bam in downloaded_bams:
        print('Processing', bam)
        process_bam(bam, out_path)
