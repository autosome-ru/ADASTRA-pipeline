import json

parameters_path = '/home/abramov/PARAMETERS/'
parallel_parameters_path = '/home/abramov/ParallelParameters/'

if __name__ == "__main__":
    with open(parameters_path + 'CELL_LINES.json', 'r') as read_file:
        d = json.loads(read_file.readline())

    keys = sorted(d.keys())
    with open(parallel_parameters_path + 'PE_parameters.cfg', 'w') as file:
        for key in keys:
            file.write(key + '\n')
