import sys
import numpy as np
import pandas as pd

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import parameters_path
from scripts.HELPERS.helpers import states


column_names = ['r', 'w', 'status', 'gof']
for BAD in states:
    ref = np.load(parameters_path + 'NBweights_ref_BAD={:.1f}.npy'.format(BAD))
    alt = np.load(parameters_path + 'NBweights_alt_BAD={:.1f}.npy'.format(BAD))
    ref_df = pd.DataFrame(columns=column_names)
    alt_df = pd.DataFrame(columns=column_names)
    for i in range(len(column_names)):
        ref_df[column_names[i]] = ref[:, i]
        alt_df[column_names[i]] = alt[:, i]
    ref_df.to_csv(parameters_path + 'weights/NBweights_ref_BAD={:.1f}.tsv'.format(BAD), sep='\t', index=False)
    alt_df.to_csv(parameters_path + 'weights/NBweights_ref_BAD={:.1f}.tsv'.format(BAD), sep='\t', index=False)
