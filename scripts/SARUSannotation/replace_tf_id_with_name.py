import os
import pandas as pd
from scripts.HELPERS.paths_for_components import results_path


def generate_uniprot_dict(uniprot_dict_file):
    u_df = pd.read_table(uniprot_dict_file)
    return pd.Series(u_df['Entry name'].values, index=u_df['Entry']).to_dict()


def main(uniprot_dict_file):
    uniprot_dict = generate_uniprot_dict(uniprot_dict_file)
    for directory in 'TF_P-values', 'TF_DICTS':
        base_dir = os.path.join(results_path, directory)
        for file in os.listdir(base_dir):
            root_ext = list(os.path.splitext(file))
            try:
                if root_ext[0].endswith('_DICT'):
                    # Fix for old format dicts
                    root_ext[0] = root_ext[0][:-5]
                new_file = uniprot_dict[root_ext[0]] + root_ext[1]
            except KeyError:
                print('No name found for given id {}'.format(root_ext[0]))
                continue
            print('Renaming {} in {}'.format(file, new_file))
            os.rename(os.path.join(base_dir, file), os.path.join(base_dir, new_file))
