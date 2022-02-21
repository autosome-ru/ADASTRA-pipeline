import os
import subprocess

from scripts.HELPERS.paths_for_components import alignments_path


dnase_bams_path = '/mnt/NAS/home/abramov/raw_alignments.GTRD/dnase'


def make_bams_list(gtrd_id):
    return [file for file in os.listdir(dnase_bams_path)
            if os.path.splitext(file)[0].split('_')[0] == gtrd_id]


def process_bam(bam, out_path):
    cmd = ['bash', 'scripts/SNPcalling.sh',
           '-Exp', os.path.join(dnase_bams_path, bam),
           '-Out', out_path]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
    if error:
        print(error)


def main(master_line):
    exp_name, align_name, read_groups, _ = master_line.split('\t')
    print('Processing', exp_name)
    out_path = os.path.join(alignments_path)
    downloaded_bams = make_bams_list(align_name)
    if downloaded_bams is None or len(downloaded_bams) < 2:
        return

    for bam in downloaded_bams:
        process_bam(bam, out_path)
