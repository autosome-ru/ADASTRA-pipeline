import json
import os.path
from scripts.HELPERS.paths_for_components import parallel_parameters_path, badmaps_dict_path

out_path = os.path.join(parallel_parameters_path, 'PE_parameters.cfg')


if __name__ == "__main__":
    with open(badmaps_dict_path, 'r') as read_file:
        d = json.loads(read_file.readline())
    keys = sorted(d.keys())
    with open(out_path, 'w') as file:
        for key in keys:

            is_empty = True
            for value in d[key]:
                if os.path.isfile(value):
                    is_empty = False
            if is_empty:
                continue
            file.write(key + '\n')
