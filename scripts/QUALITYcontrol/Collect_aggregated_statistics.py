import sys
import os
import numpy as np
import pandas as pd

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import parameters_path
from scripts.HELPERS.helpers import states

if __name__ == '__main__':
    for file_name in os.listdir(os.path.expanduser('~/DATA/TF_P-values/')):
        stats = pd.read_table('~/DATA/TF_P-values/' + file_name)
        