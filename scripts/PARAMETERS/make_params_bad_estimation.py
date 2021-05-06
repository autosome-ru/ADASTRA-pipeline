import json
import os.path

from scripts.HELPERS.paths import get_ending
from scripts.HELPERS.paths_for_components import parallel_parameters_path, badmaps_dict_path
from scripts.HELPERS.helpers import read_synonims

out_path = os.path.join(parallel_parameters_path, 'BE_parameters.cfg')


def in_cosmic(cell_line_name, cosmic_names):
    if cell_line_name == '22RV1__prostate_carcinoma_@GSE120738':
        return False
    if cell_line_name.split('@')[0] not in cosmic_names:
        return True
    return False


def main(only_cosmic=False):
    cosmic_names, _ = read_synonims()
    with open(badmaps_dict_path, 'r') as read_file:
        d = json.loads(read_file.readline())
    keys = sorted(d.keys())
    print(len(keys))
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
            if only_cosmic:
                with open(os.path.join(parallel_parameters_path, 'BE_states_parameters.cfg'), 'w') as f:
                    for caic in range(3, 6):
                        for state_sign in (
                                'int_6',
                                'full_5_and_6',
                                'full_6_but_1.33',
                                'full_6_but_2.5',
                                'full_6'
                        ):
                            file.write("{},{},{}\n".format(key, state_sign, caic))
            else:
                file.write(key + '\n')



if __name__ == '__main__':
    main()
