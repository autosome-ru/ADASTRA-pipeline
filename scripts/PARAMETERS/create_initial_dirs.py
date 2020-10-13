import os
from scripts.HELPERS.paths_for_components import badmaps_path, results_path, alignments_path
from scripts.HELPERS.paths import get_release_stats_path


def check_and_create_dir(dir_name):
    if not os.path.isdir(dir_name):
        if os.path.isfile(dir_name):
            raise AssertionError("Can't create dir {} (file with such name exists)".format(dir_name))
        os.mkdir(dir_name)


def main():
    check_and_create_dir(alignments_path)
    # Dirs for BADcalling
    check_and_create_dir(badmaps_path)
    check_and_create_dir(os.path.join(badmaps_path, 'merged_vcfs'))
    check_and_create_dir(os.path.join(badmaps_path, 'CAIC'))
    # Dirs for ASBcalling
    check_and_create_dir(results_path)
    check_and_create_dir(get_release_stats_path())
    for what_for in 'TF', 'CL':
        check_and_create_dir(os.path.join(results_path, what_for + '_DICTS'))
        check_and_create_dir(os.path.join(results_path, what_for + '_P-values'))


if __name__ == '__main__':
    main()
