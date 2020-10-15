import os
import sys

from scripts.HELPERS.paths_for_components import parallel_parameters_path, correlation_path


def main():
    with open(os.path.join(parallel_parameters_path, 'CS_parameters.cfg'), 'w') as file:
        for file_name in os.listdir(correlation_path):
            if os.path.isdir(os.path.join(correlation_path, file_name)) and file_name.endswith('_tables'):
                for file_name2 in os.listdir(os.path.join(correlation_path, file_name)):
                    file.write(file_name2 + '\n')
                break


if __name__ == "__main__":
    main()
