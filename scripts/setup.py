from setuptools import setup, find_packages

setup(
    name='ADASTRA',
    packages=find_packages(),
    version='1.0.0',
    entry_points={
        'console_scripts': [
            'ADASTRA badmaps_dict = PARAMETERS.make_badmaps_dict:main',
            'ADASTRA sort_cols = SNPcalling.sort_columns:main',
            'ADASTRA init_dirs = PARAMETERS.create_initial_dirs:main',
            'ADASTRA aggregation_dict = PARAMETERS.make_aggregation_dict:main',
            'ADASTRA make_paths = PARAMETERS.make_exp_paths_from_master_list:main',
            'ADASTRA badmaps_params = PARAMETERS.make_params_bad_estimation:main',
            'ADASTRA aggregation_params = PARAMETERS.make_params_aggregation:main',
            'ADASTRA sort_params = PARAMETERS.sort_params:main',
            'ADASTRA check_pos_peaks = PEAKannotation.check_pos_peaks:main',
            'ADASTRA annotate_peaks = PEAKannotation.annotate:main'

            'ADASTRA vcf_merge = BADcalling.VCFMerger:main',
            'ADASTRA bad_call = BADcalling.BADestimation:main',
            'ADASTRA bad_annotation = ASBcalling.BAD_annotation:main',
            'ADASTRA collect_ref_bias = FITnoise.collect_ref_bias_statistics:main',
            'ADASTRA fit_neg_bin = FITnoise.fit_negative_binom_with_weights:main',
            'ADASTRA neg_bin_p = ASBcalling.NBpcounter:main',
            'ADASTRA aggregation = ASBcalling.Aggregation:main'
        ],
    },
    author="Sergey Abramov, Alexandr Boytsov",
    install_requires=[
        'numpy>=1.18.0',
        'pandas>=1.0.4',
        'matplotlib>=3.2.1',
        'seaborn>=0.10.1',
        'docopt>=0.6.2',
        'requests>=2.24.0',
        'statsmodels>=0.11.1'
    ],
    python_requires='>=3.6',
)
