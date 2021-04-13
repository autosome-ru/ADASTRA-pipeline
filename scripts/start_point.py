"""
Usage:
            adastra badmaps_dict
            adastra sort_cols
            adastra init_dirs
            adastra aggregation_dict
            adastra make_paths --mode <mode>
            adastra badmaps_params
            adastra aggregation_params --for <for>
            adastra annotation_params
            adastra correlation_params
            adastra sort_params
            adastra check_pos_peaks --peak <path> --out <path> --type <type>
            adastra annotate_peaks --base <path>
            adastra vcf_merge --group <group>
            adastra bad_call --group <group>
            adastra bad_annotation --base <path>
            adastra collect_ref_bias [stats] [--suffix <suffix>] [--cell-type <name>]
            adastra fit_neg_bin
            adastra neg_bin_p --base <path>
            adastra aggregation --for <for> --name <name>
            adastra annotate_snps_for_correlation --base <path>
            adastra cosmic_correlation --base <path>
            adastra join_correlation_threads
            adastra collect_release_stats
            adastra weights_to_df
            adastra extract_sarus_data --name <name> --motif-len <int>
            adastra annotate_table_with_sarus --name <name> --motif-len <int>
            adastra annotate_with_phenotypes
            adastra extract_context
            adastra count_p <exp> <aligns>
            adastra -h | --help

Arguments:
    <mode>     Name of the mode [default: badmaps]
    <name>     Name of params(TF or CL)
    <type>     Peak type (gem, sissrs, peaks, cpics)
    <path>     Path to file
    <suffix>   Suffix for stats file
    <int>      Positive integer

Options:
    -h, --help                  Show help.
    --mode=<mode>               Mode for make_paths
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
    --uniprot-file=<path>       Path to file with uniprot conversion
"""
import time

from docopt import docopt
from babachi import BADEstimation

from .HELPERS.helpers import segmentation_states
from .HELPERS.paths import create_merged_vcf_path_function, create_badmaps_path_function


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
        main(args['--mode'])
    elif args['badmaps_params']:
        from .PARAMETERS.make_params_bad_estimation import main
        main()
    elif args['annotation_params']:
        from .PARAMETERS.make_params_annotation import main
        main()
    elif args['correlation_params']:
        from .PARAMETERS.make_params_correlation import main
        main()
    elif args['aggregation_params']:
        from .PARAMETERS.make_params_aggregation import main
        main(args['--for'])
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
        main(args['--group'])
    elif args['bad_call']:
        bad_group = args['--group']
        t = time.clock()
        with open(create_merged_vcf_path_function(bad_group)) as m_vcf:
            snps_collection, chromosomes_order, _ = BADEstimation.parse_input_file(m_vcf, allele_reads_tr=5)
            GS = BADEstimation.GenomeSegmentator(
                snps_collection=snps_collection,
                chromosomes_order=chromosomes_order,
                out=create_badmaps_path_function(bad_group),
                states=segmentation_states,
                b_penalty=4,
                verbose=True,
                allele_reads_tr=5,
                segmentation_mode='corrected'
            )
            GS.estimate_BAD()
        print('Total time: {} s'.format(time.clock() - t))

    elif args['bad_annotation']:
        from .ASBcalling.BAD_annotation import main
        main(args['--base'])
    elif args['collect_ref_bias']:
        from .FITnoise.collect_ref_bias_statistics import main
        if not args['stats']:
            main()
        else:
            main(args['--cell-type'], args['--suffix'], True)
    elif args['fit_neg_bin']:
        from .FITnoise.fit_negative_binom_with_weights import main
        main()
    elif args['neg_bin_p']:
        from .ASBcalling.NBpcounter import main
        main(args['--base'])
    elif args['aggregation']:
        from .ASBcalling.Aggregation import main
        main(args['--for'], args['--name'])
    elif args['annotate_snps_for_correlation']:
        from .CORRELATIONanalysis.Annotate_SNPs_with_BADmaps import main
        main(args['--base'])
    elif args['cosmic_correlation']:
        from .CORRELATIONanalysis.CorStats import main
        main(args['--base'])
    elif args['join_correlation_threads']:
        from .CORRELATIONanalysis.JoinThreads import main
        main()
    elif args['weights_to_df']:
        from .Qcontrol.neg_bin_weights_to_df import main
        main()
    elif args['collect_release_stats']:
        from .Qcontrol.check_made import main
        main()
    elif args['extract_sarus_data']:
        from .SARUSannotation.extract_sarus_data import main
        main(args['--name'], convert_motif_len_to_int(args['--motif-len']))
    elif args['annotate_table_with_sarus']:
        from .SARUSannotation.annotate_table_with_sarus import main
        main(args['--name'], convert_motif_len_to_int(args['--motif-len']))
    elif args['annotate_with_phenotypes']:
        from .PARSEphenotypes.asb_gwas_eqtl import main
        main()
    elif args['extract_context']:
        from .Qcontrol.extract_context import main
        main()
    elif args['count_p']:
        from .ASBcalling.NBpcounter import manual
        manual(args['<exp>'], args['<aligns>'])


def convert_motif_len_to_int(motif_len_string):
    if not motif_len_string:
        return None
    else:
        return int(motif_len_string)
