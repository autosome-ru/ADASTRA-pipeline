"""
Usage:
            adastra badmaps_dict
            adastra sort_cols
            adastra init_dirs
            adastra aggregation_dict
            adastra make_paths
            adastra badmaps_params
            adastra aggregation_params
            adastra sort_params
            adastra check_pos_peaks
            adastra annotate_peaks
            adastra vcf_merge
            adastra bad_call
            adastra bad_annotation
            adastra collect_ref_bias
            adastra fit_neg_bin
            adastra neg_bin_p
            adastra aggregation
            adastra -h | --help

Options:
    -h, --help                  Show help.
"""
from docopt import docopt


def main():
    args = docopt(__doc__)
    if args['badmaps_dict']:
        from scripts.PARAMETERS.make_badmaps_dict import main
        main()
    elif args['sort_cols']:
        from scripts.SNPcalling.sort_columns import main
        main()
    elif args['init_dirs']:
        from scripts.PARAMETERS.create_initial_dirs import main
        main()
    elif args['aggregation_dict']:
        from scripts.PARAMETERS.make_aggregation_dict import main
        main()
    elif args['make_paths']:
        from scripts.PARAMETERS.make_exp_paths_from_master_list import main
        main()
    elif args['badmaps_params']:
        from scripts.PARAMETERS.make_params_bad_estimation import main
        main()
    elif args['aggregation_params']:
        from scripts.PARAMETERS.make_params_aggregation import main
        main()
    elif args['sort_params']:
        from scripts.PARAMETERS.sort_params import main
        main()
    elif args['check_pos_peaks']:
        from scripts.PEAKannotation.check_pos_peaks import main
        main()
    elif args['annotate_peaks']:
        from scripts.PEAKannotation.annotate import main
        main()
    elif args['vcf_merge']:
        from scripts.BADcalling.VCFMerger import main
        main()
    elif args['bad_call']:
        from scripts.BADcalling.BADEstimation import main
        main()
    elif args['bad_annotation']:
        from scripts.ASBcalling.BAD_annotation import main
        main()
    elif args['collect_ref_bias']:
        from scripts.FITnoise.collect_ref_bias_statistics import main
        main()
    elif args['fit_neg_bin']:
        from scripts.FITnoise.fit_negative_binom_with_weights import main
        main()
    elif args['neg_bin_p']:
        from scripts.ASBcalling.NBpcounter import main
        main()
    elif args['aggregation']:
        from scripts.ASBcalling.Aggregation import main
        main()

