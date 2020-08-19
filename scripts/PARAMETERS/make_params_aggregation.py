import json
import sys
import os.path
from scripts.HELPERS.paths import open_aggregation_dict
from scripts.HELPERS.paths_for_components import parallel_parameters_path
from scripts.HELPERS.helpers import check_if_in_expected_args

out_path = os.path.join(parallel_parameters_path, 'Agr_parameters.cfg')


def main(what_for):
    check_if_in_expected_args(what_for)
    aggregation_dict_path = open_aggregation_dict(what_for)
    with open(aggregation_dict_path, 'r') as read_file:
        d = json.loads(read_file.readline())
    with open(out_path, 'w') as file:
        for key in sorted(d.keys()):
            is_empty = True
            for value in d[key]:
                if os.path.isfile(value):
                    is_empty = False
            if is_empty:
                continue
            file.write(key + '\n')


if __name__ == '__main__':
    main(sys.argv[1])
