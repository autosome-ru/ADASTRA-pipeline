import json
import os.path

from scripts.HELPERS.paths import get_ending, get_new_badmaps_dict_path
from scripts.HELPERS.paths_for_components import parallel_parameters_path, badmaps_dict_path
from scripts.HELPERS.helpers import read_synonims, default_model, get_models_list

out_path = os.path.join(parallel_parameters_path, 'BE_parameters.cfg')


def in_cosmic(cell_line_name, cosmic_names):
    # if cell_line_name == '22RV1__prostate_carcinoma_@GSE120738':
    #     return False
    if cell_line_name.split('@')[0] in cosmic_names:
        return True
    return False


def main(only_cosmic=False, remake=False):
    if remake:
        with open(get_new_badmaps_dict_path(default_model), 'r') as read_file:
            d = json.loads(read_file.readline())
    else:
        with open(badmaps_dict_path, 'r') as read_file:
            d = json.loads(read_file.readline())
    keys = sorted(d.keys())
    cosmic_names, _ = read_synonims()
    with open(out_path, 'w') as file:
        for key in keys:
            if only_cosmic and not in_cosmic(key, cosmic_names):
                continue
            is_empty = True
            for value in d[key]:
                if os.path.isfile(value + get_ending('vcf')):
                    is_empty = False
                    break
            if is_empty:
                continue
            for model in get_models_list():
                file.write("{},{}\n".format(key, model))


if __name__ == '__main__':
    main()
