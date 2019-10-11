import os

parameters_dir = '/home/abramov/ParallelParameters/'
Ploidy_path = '/home/abramov/Ploidy/'
Corelation_path = '/home/abramov/Correlation/'

file_dirs = {
    'PE_parameters.cfg': Ploidy_path,
    'ASWP_parameters.cfg': Ploidy_path,
    'CS_parameters.cfg': Corelation_path + 'Binomial_tables/',
}


def refactor_filename(line, param_file):
    if param_file == 'PE_parameters.cfg':
        return line + '.tsv'
    return line


for file_name in os.listdir(parameters_dir):
    if file_name.endswith('.cfg') and file_name != 'tmp.cfg':
        with open(parameters_dir + file_name, 'r') as file, open(parameters_dir + 'tmp.cfg', 'w') as tmp:
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
        os.remove(parameters_dir + file_name)
        print(parameters_dir + file_name)
        os.rename(parameters_dir + 'tmp.cfg', parameters_dir + file_name)