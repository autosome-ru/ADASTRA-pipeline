"""
Usage:
            adastra badmaps_dict
            adastra sort_cols
            adastra init_dirs
            adastra aggregation_dict
            adastra make_paths [--remade] --mode <mode>
            adastra badmaps_params [--remake]
            adastra aggregation_params [--remade] --for <for>
            adastra annotation_params
            adastra correlation_params [--remake]
            adastra sort_params
            adastra check_pos_peaks --peak <path> --out <path> --type <type>
            adastra annotate_peaks --base <path>
            adastra vcf_merge [--remake] --group <group>
            adastra bad_call [--remake] --group <group>
            adastra bad_annotation [--remade] --base <path>
            adastra collect_ref_bias [--remade] [stats] [--suffix <suffix>] [--cell-type <name>]
            adastra fit_neg_bin [--model <model>]
            adastra mixALime [--remade --njobs <int> --rescale-mode <rescale> --dist <dist>]
            adastra neg_bin_p --base <path>
            adastra aggregation [--remade] --for <for> --name <name>
            adastra annotate_snps_for_correlation --base <path> [--remake]
            adastra cosmic_correlation [--remake] --base <path>
            adastra join_correlation_threads [--remake]
            adastra collect_release_stats
            adastra weights_to_df
            adastra collect_roc [cell_line_wise] --group <group>
            adastra extract_sarus_data --name <name> --motif-len <int>
            adastra annotate_table_with_sarus --name <name> --motif-len <int>
            adastra annotate_with_phenotypes [--dir <path>]
            adastra extract_context
            adastra create_badmaps_filter [--njobs <int>]
            adastra apply_badmaps_filter
            adastra -h | --help

Arguments:
    <mode>     Name of the mode
    <name>     Name of params(TF or CL)
    <type>     Peak type (gem, sissrs, peaks, cpics)
    <path>     Path to file
    <suffix>   Suffix for stats file
    <int>      Positive integer
    <model>    Model to fit distribution with
    <rescale>  'none', 'single' or 'group'

Options:
    -h, --help                  Show help.
    --mode=<mode>               Mode for make_paths [default: badmaps]
    --name=<name>               Name of TF or CL
    --for=<for>                 TF or CL
    --peak=<path>               Path to peak file
    --out=<path>                Path to out file
    --type=<type>               Peak type
    --base=<path>               Path to file to annotate
    --group=<group>             Name of badmap group
    --suffix=<suffix>           Suffix for stats file
    --cell-type=<name>          Cell type name
    --motif-len=<int>           Length of the motif
    --dir=<path>                Path to directory
    --uniprot-file=<path>       Path to file with uniprot conversion
    --njobs=<int>               Number of parallel processes [default: 1]
    --rescale-mode=<rescale>    Mode of weights rescaling in mixALime
    --dist=<dist>               Dist to use. One of NB, BetaNB [default: BetaNB]
"""
import time

from docopt import docopt
from babachi import BADEstimation

from .HELPERS.helpers import get_states, create_badmaps_path_function, get_params_from_model_name, get_prior
from .HELPERS.paths import create_merged_vcf_path_function


def main():
    args = docopt(__doc__)
    if args['badmaps_dict']:
        from .PARAMETERS.make_badmaps_dict import main
        main()
    elif args['sort_cols']:
        from .SNPcalling.sort_columns import main
        main()
    elif args['init_dirs']:
        from .PARAMETERS.create_initial_dirs import main
        main()
    elif args['aggregation_dict']:
        from .PARAMETERS.make_aggregation_dict import main
        main()
    elif args['make_paths']:
        from .PARAMETERS.make_exp_paths_from_master_list import main
        main(args['--mode'], remade=args['--remade'])
    elif args['badmaps_params']:
        from .PARAMETERS.make_params_bad_estimation import main
        main(args['--remake'])
    elif args['annotation_params']:
        from .PARAMETERS.make_params_annotation import main
        main()
    elif args['correlation_params']:
        from .PARAMETERS.make_params_correlation import main
        main(args['--remake'])
    elif args['aggregation_params']:
        from .PARAMETERS.make_params_aggregation import main
        main(args['--for'], remade=args['--remade'])
    elif args['sort_params']:
        from .PARAMETERS.sort_params import main
        main()
    elif args['check_pos_peaks']:
        from .PEAKannotation.check_pos_peaks import main
        main(args['--peak'], args['--out'], args['--type'])
    elif args['annotate_peaks']:
        from .PEAKannotation.annotate import main
        main(args['--base'])
    elif args['vcf_merge']:
        from .BADcalling.VCFMerger import main
        bad_group, model = args['--group'].split(',')
        main(bad_group)
    elif args['bad_call']:
        bad_group, model = args['--group'].split(',')
        params = get_params_from_model_name(model)
        t = time.clock()
        with open(create_merged_vcf_path_function(bad_group)) as m_vcf:
            snps_collection, chromosomes_order, _ = BADEstimation.parse_input_file(m_vcf, allele_reads_tr=5)
            GS = BADEstimation.GenomeSegmentator(
                snps_collection=snps_collection,
                chromosomes_order=chromosomes_order,
                out=create_badmaps_path_function(bad_group,
                                                 valid=args['--remake'],
                                                 model=model),
                states=get_states(params['states_set']),
                b_penalty=convert_string_to_int(params['b_penalty']),
                prior=get_prior(params['states_set'], params['prior']),
                verbose=True,
                allele_reads_tr=5,
                segmentation_mode='corrected',
                atomic_region_size=600,
                chr_filter=100,
                subchr_filter=3,
                min_seg_snps=3,
                min_seg_bp=1000,
            )

            GS.estimate_BAD()
        print('Total time: {} s'.format(time.clock() - t))
    elif args['collect_roc']:
        if args['cell_line_wise']:
            from .Qcontrol.collect_cell_line_wise_data_for_ROC import main
            main(args['--group'])
        else:
            from .Qcontrol.collect_data_for_ROC import main
            main(args['--group'])
    elif args['bad_annotation']:
        from .ASBcalling.BAD_annotation import main
        main(args['--base'], remade=args['--remade'])
    elif args['collect_ref_bias']:
        from .FITnoise.collect_ref_bias_statistics import main
        if not args['stats']:
            main(remade=args['--remade'])
        else:
            main(args['--cell-type'], args['--suffix'], in_stats=True, remade=args['--remade'])
    elif args['fit_neg_bin']:
        from .FITnoise.fit_dist import main
        main()
    elif args['mixALime']:
        from .ASBcalling.calc_pval import main
        main(remade=args['--remade'],
             n_jobs=int(args['--njobs']),
             rescale_mode=args['--rescale-mode'],
             dist=args['--dist'])
    elif args['neg_bin_p']:
        from .ASBcalling.NBpcounter import main
        main(args['--base'])
    elif args['aggregation']:
        from .ASBcalling.Aggregation import main
        main(args['--for'], args['--name'], remade=args['--remade'])
    elif args['annotate_snps_for_correlation']:
        from .CORRELATIONanalysis.Annotate_SNPs_with_BADmaps import main
        main(args['--base'], remake=args['--remake'])
    elif args['cosmic_correlation']:
        from .CORRELATIONanalysis.CorStats import main
        main(args['--base'], remake=args['--remake'])
    elif args['join_correlation_threads']:
        from .CORRELATIONanalysis.JoinThreads import main
        main(args['--remake'])
    elif args['weights_to_df']:
        from .Qcontrol.neg_bin_weights_to_df import main
        main()
    elif args['collect_release_stats']:
        from .Qcontrol.check_made import main
        main()
    elif args['extract_sarus_data']:
        from .SARUSannotation.extract_sarus_data import main
        main(args['--name'], convert_string_to_int(args['--motif-len']))
    elif args['annotate_table_with_sarus']:
        from .SARUSannotation.annotate_table_with_sarus import main
        main(args['--name'], convert_string_to_int(args['--motif-len']))
    elif args['annotate_with_phenotypes']:
        from .PARSEphenotypes.asb_gwas_eqtl import main
        main(args['--dir'] if args['--dir'] else None)
    elif args['extract_context']:
        from .Qcontrol.extract_context import main
        main()
    elif args['create_badmaps_filter']:
        from scripts.BADMAPSfilter.construct_badmaps_filter import main
        main(20, 50, int(args['--njobs']))
    elif args['apply_badmaps_filter']:
        from scripts.BADMAPSfilter.apply_filter import main
        main()


def convert_string_to_int(string):
    if not string:
        return None
    else:
        return int(string)
