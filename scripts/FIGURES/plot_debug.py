import os
import sys
import json
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd

if __name__ == '__main__':
    with open(os.path.expanduser('~/TF_FC/debug_dict_fit_p.json')) as j:
        d = json.loads(j.readline())
    for tr in d:
        x = np.array(d[tr]['x'])

        fig, ax = plt.subplots(figsize=(10, 8))
        plt.scatter(-np.log10(x), d[tr]['y1'], color='blue')
        plt.scatter(-np.log10(x), d[tr]['y2'], color='red')
        plt.xlim(-1, 20)
        plt.title('CTCF min(ref,alt)_tr = {}'.format(tr))
        plt.grid(True)
        plt.xlabel('-log10 p-value')
        plt.ylabel('cumulative count')
        plt.savefig(os.path.expanduser('~/TF_FC/CTCF_count_min_tr={}').format(tr))
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(10, 8))
        plt.scatter(-np.log10(x), [b/(b+r) for b, r in zip(d[tr]['y1'], d[tr]['y2'])], color='black')
        plt.xlim(-1, 20)
        plt.ylim(0.85, 1)
        plt.title('CTCF min(ref,alt)_tr = {}'.format(tr))
        plt.grid(True)
        plt.xlabel('-log10 p-value')
        plt.ylabel('cumulative blue fraction')
        plt.savefig(os.path.expanduser('~/TF_FC/CTCF_fraction_min_tr={}').format(tr))
        plt.close(fig)
