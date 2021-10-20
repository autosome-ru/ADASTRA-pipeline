import os
import sys

from scripts.HELPERS.paths import get_correlation_path
from scripts.HELPERS.paths_for_components import parallel_parameters_path


def main(remake=False):
    with open(os.path.join(parallel_parameters_path, 'CS_parameters.cfg'), 'w') as file:
        for file_name in os.listdir(get_correlation_path()):
            if os.path.isdir(os.path.join(get_correlation_path(), file_name)) and file_name.endswith('_tables{}'.format('_filtered' if remake else '')):
                for file_name2 in os.listdir(os.path.join(get_correlation_path(), file_name)):
                    file.write(file_name2 + '\n')
                break


if __name__ == "__main__":
    main()
