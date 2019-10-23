import os

parallel_parameters_dir = '/home/abramov/ParallelParameters/stats/'
Ploidy_path = '/home/abramov/PloidyForRelease/'
Correlation_path = '/home/abramov/Correlation/'

file_dirs = {
    'PE_parameters.cfg': Ploidy_path,
    'ASWP_parameters.cfg': Ploidy_path,
    'CS_parameters.cfg': Correlation_path + 'Binomial_tables/',
}


def refactor_filename(line, param_file):
    if param_file == 'PE_parameters.cfg':
        return line + '.tsv'
    return line


for file_name in os.listdir(parallel_parameters_dir):
    if file_name in file_dirs:
        with open(parallel_parameters_dir + file_name, 'r') as file, open(parallel_parameters_dir
                                                                          + 'tmp.cfg', 'w') as tmp:
            assert file_name != 'tmp.cfg'
            lines = []
            for line in file:
                lines.append(line.strip())
            f = lambda x: refactor_filename(file_dirs[file_name] + x,
                            param_file=file_name)
            lines = sorted(lines,
                            key=lambda x: os.path.getsize(f(x)) if os.path.isfile(f(x)) else -1,
                            reverse=True)
            for line in lines:
                tmp.write(line + '\n')
        os.remove(parallel_parameters_dir + file_name)
        print(parallel_parameters_dir + file_name)
        os.rename(parallel_parameters_dir + 'tmp.cfg', parallel_parameters_dir + file_name)