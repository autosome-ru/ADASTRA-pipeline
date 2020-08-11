import os
from scripts.HELPERS.paths_for_components import parallel_parameters_path, badmaps_path

file_dirs = {
    'BE_parameters.cfg': badmaps_path,
}


def refactor_filename(line, param_file):
    if param_file == 'BE_parameters.cfg':
        return line + '.tsv'
    return line


for file_name in os.listdir(parallel_parameters_path):
    if file_name in file_dirs:
        with open(os.path.join(parallel_parameters_path + file_name), 'r') as file, open(
                os.path.join(parallel_parameters_path, 'tmp.cfg'), 'w') as tmp:
            assert file_name != 'tmp.cfg'
            lines = []
            for line in file:
                lines.append(line.strip())
            f = lambda x: refactor_filename(os.path.join(badmaps_path, x),
                                            param_file=file_name)
            lines = sorted(lines,
                           key=lambda x: os.path.getsize(f(x)) if os.path.isfile(f(x)) else -1,
                           reverse=True)
            for line in lines:
                tmp.write(line + '\n')
        os.remove(os.path.join(parallel_parameters_path, file_name))
        print("Made parameters {}".format(file_name))
        os.rename(os.path.join(parallel_parameters_path + 'tmp.cfg'), os.path.join(parallel_parameters_path, file_name))
