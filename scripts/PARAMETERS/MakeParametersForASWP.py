import os

Ploidy_path = '/home/abramov/Ploidy/'
parallel_parameters_path = '/home/abramov/ParallelParameters/stats/'

if __name__ == "__main__":
    with open(parallel_parameters_path + 'ASWP_parameters.cfg', 'w') as file:
        for file_name in os.listdir(Ploidy_path):
            if not os.path.isfile(Ploidy_path + file_name):
                continue
            file.write(file_name + '\n')
