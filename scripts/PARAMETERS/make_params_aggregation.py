import json
import sys
import os.path

from scripts.HELPERS.paths import get_aggregation_dict_path, get_merged_badmaps_dict_path
from scripts.HELPERS.paths_for_components import parallel_parameters_path
from scripts.HELPERS.helpers import check_if_in_expected_args, split_ext_recursive, make_reverse_dict, is_valid

out_path = os.path.join(parallel_parameters_path, 'Agr_parameters.cfg')


def main(what_for, remade=True):
    check_if_in_expected_args(what_for)
    aggregation_dict_path = get_aggregation_dict_path(what_for)
    with open(get_merged_badmaps_dict_path(remade=remade), "r") as read_file:
        d = json.load(read_file)
        rev_d = make_reverse_dict(d)
    with open(aggregation_dict_path, 'r') as read_file:
        d = json.loads(read_file.readline())
    with open(out_path, 'w') as file:
        for key in sorted(d.keys()):
            is_empty = True
            for value in d[key]:
                if os.path.isfile(value) and is_valid(split_ext_recursive(value), rev_d, remade=remade):
                    is_empty = False
            if is_empty:
                continue
            file.write(key + '\n')


if __name__ == '__main__':
    main(sys.argv[1])
