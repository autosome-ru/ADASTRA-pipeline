import os
import sys

from scripts.HELPERS.paths_for_components import parallel_parameters_path, badmaps_path


def main():
    with open(os.path.join(parallel_parameters_path, 'ASWP_parameters.cfg'), 'w') as file:
        for file_name in os.listdir(os.path.join(badmaps_path, 'merged_vcfs')):
            badmap = os.path.join(badmaps_path, 'merged_vcfs/', file_name)
            if not os.path.isfile(badmap):
                continue
            if not badmap.endswith('.tsv'):
                continue
            file.write(file_name + '\n')


if __name__ == "__main__":
    main()
