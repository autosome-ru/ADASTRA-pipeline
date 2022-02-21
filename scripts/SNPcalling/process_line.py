import os
import subprocess

from scripts.HELPERS.paths_for_components import alignments_path


dnase_bams_path = ''

def make_bams_map():
    result = {}
    for file in os.listdir(dnase_bams_path):
        file_name, _ = os.path.splitext(file)[0].split('_')
        result.setdefault(file_name, []).append(file)
    return result

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
    out_path = os.path.join(alignments_path)
    downloaded_bams = make_bams_map().get(exp_name, None)
    if downloaded_bams is None or len(downloaded_bams) == 1:
        return
    for bam in downloaded_bams:
        process_bam(bam, out_path)
