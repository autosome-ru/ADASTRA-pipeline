import sys
import numpy as np
import os
import json
from scipy import stats
import pandas as pd
sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.FIGURES import style_config
from matplotlib import pyplot as plt
import seaborn as sns

if __name__ == '__main__':
    with open(os.path.expanduser('~/phantom_style/overall_statistics.json')) as j:
        d = json.loads(j.readline())

    print(d.keys())

    sns.set(font_scale=1.4, style="ticks", font="lato", palette=('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2',
                                                             '#D55E00', '#CC79A7'))
    sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})
    plt.rcParams['font.weight'] = "medium"
    plt.rcParams['axes.labelweight'] = 'medium'
    plt.rcParams['figure.titleweight'] = 'medium'
    plt.rcParams['axes.titleweight'] = 'medium'
    plt.rcParams['figure.figsize'] = 6, 4
    plt.rcParams["legend.framealpha"] = 1
    plt.rcParams['axes.xmargin'] = 0
    plt.rcParams['axes.ymargin'] = 0
    plt.rcParams["legend.framealpha"] = 1


    #Total SNP calls
    vals = list(d['SNP_calls']['TF'].values())
    vals = sorted(vals, reverse=True)
    fig, ax = plt.subplots()
    fig.tight_layout(pad=1.5)
    # ax.margins(x=0, y=0)
    plt.stackplot(range(len(vals)), vals, alpha=0.7, linewidth=0.7, color='C4')
    # plt.yticks(range(1, 8), [''])
    plt.ylabel('Number of SNP calls')
    plt.yscale('log')
    plt.xlabel('TFs sorted by number of datasets')
    plt.savefig(os.path.expanduser('~/phantom_style/Figure_AS_1_total_snp_barplot_TF.svg'), dpi=300)
    # plt.show()
    plt.close(fig)

    vals = list(d['SNP_calls']['CL'].values())
    vals = sorted(vals, reverse=True)
    fig, ax = plt.subplots()
    fig.tight_layout(pad=1.5)
    # ax.margins(x=0, y=0)
    plt.stackplot(range(len(vals)), vals, alpha=0.7, linewidth=0.7, color='C4')
    # plt.yticks(range(1, 8), [''])

    plt.ylabel('Number of SNP calls')
    plt.yscale('log')
    plt.xlabel('Cell types sorted by number of datasets')
    plt.savefig(os.path.expanduser('~/phantom_style/Figure_AS_1_total_snp_barplot_CL.svg'), dpi=300)
    # plt.show()
    plt.close(fig)

    #Datasets
    vals = list(d['datasets']['TF'].values())
    vals = sorted(vals, reverse=True)
    fig, ax = plt.subplots()
    fig.tight_layout(pad=1.5)
    # ax.margins(x=0, y=0)
    plt.stackplot(range(len(vals)), vals, alpha=0.7, linewidth=0.7, color='C4')
    # plt.yticks(range(1, 8), [''])
    plt.ylabel('Number of datasets')
    plt.yscale('log')
    plt.ylim(0.7, 405)
    plt.xlabel('TFs sorted by number of datasets')
    plt.savefig(os.path.expanduser('~/phantom_style/Figure_AS_1_datasets_barplot_TF.svg'), dpi=300)
    # plt.show()
    plt.close(fig)

    vals = list(d['datasets']['CL'].values())
    vals = sorted(vals, reverse=True)
    fig, ax = plt.subplots()
    fig.tight_layout(pad=1.5)
    # ax.margins(x=0, y=0)
    plt.stackplot(range(len(vals)), vals, alpha=0.7, linewidth=0.7, color='C4')
    # plt.yticks(range(1, 8), [''])
    plt.ylabel('Number of datasets')
    plt.yscale('log')
    plt.ylim(0.7, 405)
    plt.xlabel('Cell types sorted by number of SNP calls')
    plt.savefig(os.path.expanduser('~/phantom_style/Figure_AS_1_datasets_barplot_CL.svg'), dpi=300)
    # plt.show()
    plt.close(fig)

    #Unique SNPS
    vals = []
    vals_asb = []
    for key in d['unique_SNPs']['TF']:
        vals.append(d['unique_SNPs']['TF'][key] - d['unique_asb']['TF'][key])
        vals_asb.append(d['unique_asb']['TF'][key])
    vals, vals_asb = list(zip(*sorted(list(zip(vals, vals_asb)), key=lambda x: x[0], reverse=True)))
    fig, ax = plt.subplots()
    fig.tight_layout(pad=1.5)
    # ax.margins(x=0, y=0)
    plt.stackplot(range(len(vals)), vals, alpha=0.7, linewidth=0.7, color='#CCCCCC', labels=['non-ASB'])
    plt.stackplot(range(len(vals)), vals_asb, alpha=0.7, linewidth=0.7, color='C0', labels=['ASB'])
    # plt.yticks(range(1, 8), [''])
    plt.ylabel('Number of rsSNPs')
    plt.yscale('log')
    plt.xlabel('TFs sorted by number of rsSNPs')
    plt.legend(loc='upper right')
    plt.savefig(os.path.expanduser('~/phantom_style/Figure_AS_1_unique_snp_barplot_TF.svg'), dpi=300)
    # plt.show()

    vals = []
    vals_asb = []
    for key in d['unique_SNPs']['CL']:
        vals.append(d['unique_SNPs']['CL'][key] - d['unique_asb']['CL'][key])
        vals_asb.append(d['unique_asb']['CL'][key])
    vals, vals_asb = list(zip(*sorted(list(zip(vals, vals_asb)), key=lambda x: x[0], reverse=True)))
    fig, ax = plt.subplots()
    fig.tight_layout(pad=1.5)
    # ax.margins(x=0, y=0)
    plt.stackplot(range(len(vals)), vals, alpha=0.7, linewidth=0.7, color='#CCCCCC', labels=['non-ASB'])
    plt.stackplot(range(len(vals)), vals_asb, alpha=0.7, linewidth=0.7, color='C0', labels=['ASB'])
    # plt.yticks(range(1, 8), [''])
    plt.ylabel('Number of rsSNPs')
    plt.yscale('log')
    plt.xlabel('Cell types sorted by number of rsSNPs')
    plt.legend(loc='upper right')
    plt.savefig(os.path.expanduser('~/phantom_style/Figure_AS_1_unique_snp_barplot_CL.svg'), dpi=300)
    # plt.show()
