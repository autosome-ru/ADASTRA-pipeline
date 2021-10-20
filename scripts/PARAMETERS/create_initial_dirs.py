import os

from scripts.HELPERS.helpers import check_and_create_dir
from scripts.HELPERS.paths_for_components import badmaps_path, results_path, alignments_path
from scripts.HELPERS.paths import get_release_stats_path, get_correlation_path, get_heatmap_data_path, \
    get_badmaps_path_by_validity, stage_dict, get_dir_by_stage


def main():
    check_and_create_dir(alignments_path)

    # Dirs for ASBcalling
    check_and_create_dir(results_path)
    for stage in stage_dict:
        check_and_create_dir(get_dir_by_stage(stage))
    check_and_create_dir(get_release_stats_path())
    for what_for in 'TF', 'CL':
        check_and_create_dir(os.path.join(results_path, what_for + '_DICTS'))
        check_and_create_dir(os.path.join(results_path, what_for + '_P-values'))
    # Dirs for BADcalling
    check_and_create_dir(badmaps_path)
    check_and_create_dir(os.path.join(badmaps_path, 'merged_vcfs'))
    check_and_create_dir(get_correlation_path())
    check_and_create_dir(get_heatmap_data_path())
    for validity in True, False:
        check_and_create_dir(get_badmaps_path_by_validity(validity))
        check_and_create_dir(os.path.join(get_badmaps_path_by_validity(validity), 'CAIC'))


if __name__ == '__main__':
    main()
