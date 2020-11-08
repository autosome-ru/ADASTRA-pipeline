import os


def main(directory, download_path):
    download_path_list = download_path.split(' ')
    with open(os.path.join(directory, 'bam_list.txt'), 'w') as out:
        for file in download_path_list:
            out.write(os.path.join(directory, os.path.basename(file)) + '\n')
