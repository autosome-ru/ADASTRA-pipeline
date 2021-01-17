import pandas as pd
import os
import scipy.stats as st
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns
import numpy as np

from scripts.HELPERS.helpers import get_states
states = get_states('all_5')


def get_p_value(row):
    print(row['seg_id'])
    p = 1/(row['BAD'] + 1)
    n = row['ref'] + row['alt']
    dist = st.binom(n=n, p=p)
    cdf = dist.cdf
    x = min(row['ref'], row['alt'])
    return (cdf(x) - cdf(4) + cdf(n-5) - cdf(n-x-1)) / (cdf(n-5) - cdf(4)) if x < n/2 else 1


def extract_biosamples(name):
    if name.find('ENC') != -1:
        return name[name.find('ENC'):-5]
    if name.find('GSE') != -1:
        return name[name.find('GSE'):-4]


def extract_cells(name):
    name = name.split('/')[-1]
    for pref in ['K562', 'MCF7', 'NCI-H128', 'cranial_neural_crest_cells']:
        if name.startswith(pref):
            return pref


get_BAD_label = {
    1: '1',
    4/3: '4/3',
    3/2: '3/2',
    2: '2',
    5/2: '5/2',
    3: '3',
    4: '4',
    5: '5',
    6: '6',
}


files = [os.path.expanduser('~/Desktop/snps/{}'.format(x)) for x in
                            [
                                'K562__myelogenous_leukemia__-labs-michael-snyder---biosamples-ENCBS389PVA-.tsv',
                                'K562__myelogenous_leukemia__GSE81009.tsv',
                                'MCF7__Invasive_ductal_breast_carcinoma__-labs-michael-snyder---biosamples-ENCBS773JGQ-.tsv',
                                'NCI-H128__small_cell_lung_carcinoma__GSE41105.tsv',
                                'MCF7__Invasive_ductal_breast_carcinoma__GSE78113.tsv',
                                'cranial_neural_crest_cells_GSE70751.tsv',
                             ]
         ]

for ind, file in enumerate(files):
    with open(file) as f:
        line = f.readline()
        print(line)
        aligns = line.strip().split('@')[-1].split(',')

    print(file)
    print(aligns)

    sns.set(font_scale=1, style="ticks", font="lato")
    sns.set_palette((
        '#56B4E9',
        '#0072B2',
        '#009E73',
        '#E69F00',
        '#F0E442',
        '#D55E00',
        '#505050',
        '#CC79A7'))
    sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})
    plt.rcParams['font.weight'] = "medium"
    plt.rcParams['axes.labelweight'] = 'medium'
    plt.rcParams['figure.titleweight'] = 'medium'
    plt.rcParams['axes.titleweight'] = 'medium'
    # plt.rcParams['figure.figsize'] = 6, 3 * len(aligns)
    plt.rcParams['figure.figsize'] = 8, 8
    plt.rcParams["legend.framealpha"] = 1
    plt.rcParams["legend.framealpha"] = 1

    if os.path. isfile(file + '.pv'):
        t = pd.read_table(file + '.pv')
    else:
        t = pd.read_table(file, header=None, comment='#')
        t.columns = ['chr', 'pos', 'ref', 'alt', 'BAD'] + ['Q{:.2f}'.format(BAD) for BAD in states] + ['snps_n', 'sumcov', 'dataset', 'seg_id']
        # t = t[t['chr'] == 'chr1']
        # t = t[t['ref'] + t['alt'] >= 15]
        t['p_value'] = t.apply(get_p_value, axis=1)
        t.to_csv(file + '.pv', sep='\t', index=False)
    print('q', np.quantile(t['p_value'], 0.05))
    q_seg_wise = [np.quantile(t[t['seg_id'] == i]['p_value'], 0.05) for i in list(set(t['seg_id']))]
    print('q_list', q_seg_wise)
    q_data_wise = [np.quantile(t[t['dataset'] == i]['p_value'], 0.05) for i in list(set(t['dataset']))]
    print('q_list_datas', q_data_wise)

    # for dataset in list(set(t['dataset'])):
    #     values = [np.quantile(t[t['dataset'] == dataset]['p_value'], proc/100) for proc in range(1, 51)]
    #     plt.scatter(x=[proc/100 for proc in range(1, 51)], y=values, color='C{}'.format(ind+1))

    cov_filter = 10
    t = t[t['ref'] + t['alt'] >= cov_filter]

    xs = [proc/1000 for proc in range(1, 401)]
    values = [np.quantile(t['p_value'], x) for x in xs]
    plt.plot(values, [xs[i] for i in range(len(xs))], color='C{}'.format(ind+1), label='{}@{}'.format(extract_cells(file), extract_biosamples(file)))

plt.plot([0, 0.4], [0, 0.4], label='uniform')
plt.title('Empirical P-value CDFs\nCoverage filter: >={}'.format(cov_filter))
plt.xlabel('P-value')
plt.ylabel('CDF')
plt.grid()
plt.legend()
# plt.savefig(os.path.expanduser('~/Desktop/segval/CDF@{}.png'.format(cov_filter)))
plt.show()


    # plt.hist(q_data_wise, range=(0, 1), bins=50)
    # plt.title('{}@{}'.format(extract_cells(file), extract_biosamples(file)))
    # plt.show()
    #
    # if not os.path.isdir(os.path.expanduser('~/Desktop/segval/{}@{}'.format(extract_cells(file), extract_biosamples(file)))):
    #     os.mkdir(os.path.expanduser('~/Desktop/segval/{}@{}'.format(extract_cells(file), extract_biosamples(file))))
    #
    # for segment_id in list(set(t['seg_id'])):
    #     fig, (*axs, ) = plt.subplots(len(aligns), 1)
    #     for i, dataset in enumerate(aligns):
    #         q = t[(t['dataset'] == dataset) & (t['seg_id'] == segment_id)]
    #         sc = sum(q['ref'] + q['alt'])
    #         axs[i].hist(q['p_value'], range=(0, 1), bins=20)
    #         axs[i].set_title('{}, sumcov: {}, avgcov: {:.0f}'.format(dataset, sc, sc / len(q.index) if len(q.index) > 0 else 0))
    #     plt.suptitle('{}@{}\nSeg {}, BAD: {}, snps: {}'.format(extract_cells(file), extract_biosamples(file), segment_id, get_BAD_label[t[t['seg_id'] == segment_id]['BAD'].unique()[0]], len(t[t['seg_id'] == segment_id].index)), fontsize=18)
    #     plt.savefig(os.path.expanduser('~/Desktop/segval/{}@{}/seg_{}.png'.format(extract_cells(file), extract_biosamples(file), segment_id)))
    #     plt.close(fig)
    #
    # fig, (*axs,) = plt.subplots(len(aligns), 1)
    # for i, dataset in enumerate(aligns):
    #     q = t[t['dataset'] == dataset]
    #     sc = sum(q['ref'] + q['alt'])
    #     axs[i].hist(q['p_value'], range=(0, 1), bins=20)
    #     axs[i].set_title(
    #         '{}, sumcov: {}, avgcov: {:.0f}'.format(dataset, sc, sc / len(q.index) if len(q.index) > 0 else 0))
    # plt.suptitle('{}@{}\nAll segments, snps: {}'.format(extract_cells(file), extract_biosamples(file),
    #                                                        len(t.index)), fontsize=18)
    # plt.savefig(os.path.expanduser(
    #     '~/Desktop/segval/{}@{}/seg_all.png'.format(extract_cells(file), extract_biosamples(file))))
    # plt.close(fig)
