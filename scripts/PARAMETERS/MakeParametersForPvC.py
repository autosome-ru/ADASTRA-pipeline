import json
import os.path

parameters_path = '/home/abramov/PARAMETERS/'
parallel_parameters_path = '/home/abramov/ParallelParameters/'

if __name__ == "__main__":
    with open(parameters_path + 'REVERSE_CELL_LINES.json', 'r') as read_file:
        d = json.loads(read_file.readline())

    keys = sorted(d.keys())
    with open(parallel_parameters_path + 'PvC_parameters.cfg', 'w') as file:
        for key in keys:
            full_path = ".".join(key.split(".")[:-2])
            if os.path.isfile(key) and os.path.isfile(full_path + "_table_annotated.txt"):
                file.write(full_path + '\n')
