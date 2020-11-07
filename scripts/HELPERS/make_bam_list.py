import sys
import os

if __name__ == '__main__':
    directory = sys.argv[1]
    download_path = sys.argv[2]
    download_path_list = download_path.split(' ')
    with open(os.path.join(directory, 'bam_list.txt'), 'w') as out:
        for file in download_path_list:
            out.write(os.path.join(directory, os.path.basename(file)) + '\n')
