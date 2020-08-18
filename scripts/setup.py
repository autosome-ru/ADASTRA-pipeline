from setuptools import setup, find_packages

setup(
    name='ADASTRA',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'badmaps_dict = PARAMETERS.make_badmaps_dict:main',
            'sort_cols = SNPcalling.sort_columns:main',
            'init_dirs = PARAMETERS.create_initial_dirs:main',
            'aggregation_dict = PARAMETERS.make_aggregation_dict:main',
            'make_paths = PARAMETERS.make_exp_paths_from_master_list:main',
            'badmaps_params = PARAMETERS.make_params_bad_estimation:main',
            'aggregation_params = PARAMETERS.make_params_aggregation:main',
            'sort_params = PARAMETERS.sort_params:main',
            'check_pos_peaks = PEAKannotation.check_pos_peaks:main',
            'annotate_peaks = PEAKannotation.annotate:main'

            'vcf_merge = BADcalling.VCFMerger:main',
            'bad_call = BADcalling.BADestimation:main',
            'bad_annotation = ASBcalling.BAD_annotation:main',
            'collect_ref_bias = FITnoise.collect_ref_bias_statistics:main',
            'fit_neg_bin = FITnoise.fit_negative_binom_with_weights:main',
            'neg_bin_p = ASBcalling.NBpcounter:main',
            'aggregation = ASBcalling.Aggregation:main'
        ],
    },
    author="Sergey Abramov, Alexandr Boytsov",
    install_requires=[
        'numpy>=1.18.0',
        'pandas>=1.0.4',
        'matplotlib>=3.2.1',
        'seaborn>=0.10.1',
        'docopt>=0.6.2'
    ],
    python_requires='>=3.6',
)