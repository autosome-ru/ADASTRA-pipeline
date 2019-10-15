import json
import os.path

if __name__ == "__main__":
    JSON_path = '/home/abramov/PLOIDYcalling/REVERSE_CELL_LINES.json'
    out_path = '/home/abramov/ParallelParameters/PvC_parameters.cfg'

    with open(JSON_path, 'r') as read_file:
        d = json.loads(read_file.readline())
    keys = sorted(d.keys())

    with open(out_path, 'w') as file:
        for key in keys:
            full_path = ".".join(key.split(".")[:-2])
            if os.path.isfile(key) and os.path.isfile(full_path + "_table_annotated.txt"):
                file.write(full_path + '\n')
