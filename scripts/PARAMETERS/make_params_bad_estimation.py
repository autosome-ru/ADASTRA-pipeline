import json
import os.path

from scripts.HELPERS.paths import get_ending, get_new_badmaps_dict_path
from scripts.HELPERS.paths_for_components import parallel_parameters_path, badmaps_dict_path

out_path = os.path.join(parallel_parameters_path, 'BE_parameters.cfg')


def main(remake=False):
    if remake:
        with open(get_new_badmaps_dict_path(), 'r') as read_file:
            d = json.loads(read_file.readline())
    else:
        with open(badmaps_dict_path, 'r') as read_file:
            d = json.loads(read_file.readline())
    keys = sorted(d.keys())
    with open(out_path, 'w') as file:
        for key in keys:
            is_empty = True
            for value in d[key]:
                if os.path.isfile(value + get_ending('vcf')):
                    is_empty = False
            if is_empty:
                continue
            file.write(key + '\n')


if __name__ == '__main__':
    main()
