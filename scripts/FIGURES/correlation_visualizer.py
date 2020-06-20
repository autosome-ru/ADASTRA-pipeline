import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import seaborn as sns
import os

if __name__ == '__main__':
    corstats_path = '/home/sashok/cor/cor_stats_test.tsv'
    heatmapdata_path = '/home/sashok/Documents/ASB/Correlation/HeatmapData/'

    out_folder = '/home/sashok/Documents/ASB/corpics/'
    cor_folder = 'cor_total/'
    seg_folder = 'seg_total/'
    cor_seg_folder = 'cor_seg/'
    heat_folder = 'heat/'

    table_q_path = 'model_comparison_'

    cor_table = pd.read_csv(corstats_path, sep='\t')

    # cell_lines = ['MCF7_Invasive_ductal_breast_carcinoma']
    cell_lines = [
         'K562__myelogenous_leukemia_',
         'HCT-116__colon_carcinoma_',
         'MCF7__Invasive_ductal_breast_carcinoma_',
        # 'PC3_prostate_carcinoma',
        # 'LoVo_colorectal_adenocarcinoma',
        # 'HeLa_S3_cervical_adenocarcinoma',
        # 'all_lines'
    ]
    # cell_lines = ['HCT-116_colon_carcinoma']
    models = ['CAIC',
              # 'SQRT',
              ]

    # plots = ['# of segments', 'cor', 'cor-seg', 'cosm_heatmap']
    # plots = ['cor-seg', 'cor', '# of segments']
    plots = ['cor']
    # plots = ['cosm_heatmap']

    heatmap_reg = '3'

    v1 = 25000
    v2 = 150000

    v1c = 300000
    v2c = 3000000

    regs = dict(zip(['1', '2', '3'], [(0, v1), (v1, v2), (v2, None)]))
    regnames = dict(zip(['1', '2', '3'], ['snp_0-25k', 'snp25-150k', 'snp150k+']))

    short_name = dict(zip(
        ['K562__myelogenous_leukemia_', 'MCF7__Invasive_ductal_breast_carcinoma_', 'HCT-116__colon_carcinoma_',
         'PC3_prostate_carcinoma', 'LoVo_colorectal_adenocarcinoma', 'all_lines'],
        ['K562', 'MCF7', 'HCT116', 'PC3', 'LoVo', 'All lines']))

    for cell_line in cell_lines:
        table = cor_table[cor_table['#cell_line'] == cell_line]
        #table = cor_table[(cor_table['#cell_line'] != 'K562_myelogenous_leukemia') &
         #                 (cor_table['#cell_line'] !='MCF7_Invasive_ductal_breast_carcinoma')]
        print(table.info())

        with open(table_q_path + short_name[cell_line] + '.tsv', 'w') as otab:
            header = ['cell_line', 'model', 'states', 'N_penalty', 'distance_penalty', 'cor>0, snp 0-25k',
                      'cor>0, snp 25-150k', 'cor>0, snp 150k+', 'cor>0.3, snp 0-25k', 'cor>0.3, snp 25-150k',
                      'cor>0.3, snp 150k+']
            otab.write('\t'.join(header) + '\n')

            for model in models:
                print(model)
                if '# of segments' in plots:
                    fig = plt.figure(figsize=(10, 8))
                    sns.scatterplot(data=table,
                                    x='total_snps',
                                    y='number_of_segments_' + model,
                                    color='indianred')
                    plt.title('\n'.join((cell_line, model, '# of segments', 'NEW N')))
                    plt.grid(True)
                    plt.hlines(y=np.mean(table['total_regions']), xmin=0, xmax=max(table['total_snps']))
                    plt.savefig(out_folder + seg_folder + '_'.join([cell_line, model.replace('.', ','), 'segments']))

                if 'cor' in plots:
                    fig = plt.figure(figsize=(20, 16))

                    ax = plt.subplot(221)
                    sns.scatterplot(data=table,
                                    x='total_snps',
                                    y='cor_by_snp_' + model, ax=ax)
                    vals = [cell_line] + model.split('-')
                    plt.title('\n'.join((cell_line, model, 'Correlation vs # of snps')))
                    plt.grid(True)

                    plt.vlines(x=v1, ymin=-1, ymax=1)
                    plt.vlines(x=v2, ymin=-1, ymax=1)

                    for y_cut in (0, 0.3, 0.8):
                        plt.hlines(y=y_cut, xmin=0, xmax=max(table['total_snps'], default=10000))

                        value = len(table[(table['cor_by_snp_' + model] >= y_cut) & (table['total_snps'] <= v1)].index)
                        vals.append(value)
                        plt.text(x=v1 / 2, y=y_cut, va='bottom', ha='center', s=str(value))

                        value = len(table[(table['cor_by_snp_' + model] >= y_cut) & (table['total_snps'] > v1) & (
                                table['total_snps'] <= v2)].index)
                        vals.append(value)
                        plt.text(x=(v1 + v2) / 2, y=y_cut, va='bottom', ha='center', s=str(value))

                        value = len(table[(table['cor_by_snp_' + model] >= y_cut) & (table['total_snps'] > v2)].index)
                        vals.append(value)
                        plt.text(x=min(max(table['total_snps'] + v2, default=10000) / 2, max(table['total_snps'] * 0.95, default=10000)), y=y_cut,
                                 va='bottom', ha='center', s=str(value))

                    ax = plt.subplot(223)
                    sns.scatterplot(data=table,
                                    x='sum_cov',
                                    y='cor_by_snp_' + model, ax=ax)
                    vals = [cell_line] + model.split('-')
                    plt.title('\n'.join((cell_line, model, 'Correlation vs sum cover')))
                    plt.grid(True)

                    plt.vlines(x=v1c, ymin=-1, ymax=1)
                    plt.vlines(x=v2c, ymin=-1, ymax=1)

                    for y_cut in (0, 0.3, 0.8):
                        plt.hlines(y=y_cut, xmin=0, xmax=max(table['sum_cov'], default=10000))

                        value = len(table[(table['cor_by_snp_' + model] >= y_cut) & (table['sum_cov'] <= v1c)].index)
                        vals.append(value)
                        plt.text(x=v1c / 2, y=y_cut, va='bottom', ha='center', s=str(value))

                        value = len(table[(table['cor_by_snp_' + model] >= y_cut) & (table['sum_cov'] > v1c) & (
                                table['sum_cov'] <= v2c)].index)
                        vals.append(value)
                        plt.text(x=(v1c + v2c) / 2, y=y_cut, va='bottom', ha='center', s=str(value))

                        value = len(table[(table['cor_by_snp_' + model] >= y_cut) & (table['sum_cov'] > v2c)].index)
                        vals.append(value)
                        plt.text(x=min(max(table['sum_cov'] + v2c, default=10000) / 2, max(table['sum_cov'] * 0.95, default=10000)), y=y_cut,
                                 va='bottom', ha='center', s=str(value))

                    ax = plt.subplot(222)
                    sns.scatterplot(x=np.log10(1 + table['total_snps']),
                                    y=table['cor_by_snp_' + model], ax=ax)
                    vals = [cell_line] + model.split('-')
                    plt.title('\n'.join((cell_line, model, 'Correlation vs # of snps')))
                    plt.grid(True)

                    plt.vlines(x=np.log10(v1), ymin=-1, ymax=1)
                    plt.vlines(x=np.log10(v2), ymin=-1, ymax=1)

                    for y_cut in (0, 0.3, 0.8):
                        plt.hlines(y=y_cut, xmin=2, xmax=np.log10(1 + max(table['total_snps'], default=10000)))

                        value = len(table[(table['cor_by_snp_' + model] >= y_cut) & (table['total_snps'] <= v1)].index)
                        vals.append(value)
                        plt.text(x=v1 / 2, y=y_cut, va='bottom', ha='center', s=str(value))

                        value = len(table[(table['cor_by_snp_' + model] >= y_cut) & (table['total_snps'] > v1) & (
                                table['total_snps'] <= v2)].index)
                        vals.append(value)
                        plt.text(x=(v1 + v2) / 2, y=y_cut, va='bottom', ha='center', s=str(value))

                        value = len(table[(table['cor_by_snp_' + model] >= y_cut) & (table['total_snps'] > v2)].index)
                        vals.append(value)
                        plt.text(x=min(max(table['total_snps'] + v2, default=10000) / 2, max(table['total_snps'] * 0.95, default=10000)), y=y_cut,
                                 va='bottom', ha='center', s=str(value))

                    ax = plt.subplot(224)
                    sns.scatterplot(x=np.log10(1 + table['sum_cov']),
                                    y=table['cor_by_snp_' + model], ax=ax)
                    vals = [cell_line] + model.split('-')
                    plt.title('\n'.join((cell_line, model, 'Correlation vs sum cover')))
                    plt.grid(True)

                    plt.vlines(x=np.log10(v1c), ymin=-1, ymax=1)
                    plt.vlines(x=np.log10(v2c), ymin=-1, ymax=1)

                    for y_cut in (0, 0.3, 0.8):
                        plt.hlines(y=y_cut, xmin=3, xmax=np.log10(1 + max(table['sum_cov'], default=10000)))

                        value = len(table[(table['cor_by_snp_' + model] >= y_cut) & (table['sum_cov'] <= v1c)].index)
                        vals.append(value)
                        plt.text(x=v1c / 2, y=y_cut, va='bottom', ha='center', s=str(value))

                        value = len(table[(table['cor_by_snp_' + model] >= y_cut) & (table['sum_cov'] > v1c) & (
                                table['sum_cov'] <= v2c)].index)
                        vals.append(value)
                        plt.text(x=(v1c + v2c) / 2, y=y_cut, va='bottom', ha='center', s=str(value))

                        value = len(table[(table['cor_by_snp_' + model] >= y_cut) & (table['sum_cov'] > v2c)].index)
                        vals.append(value)
                        plt.text(x=min(max(table['sum_cov'] + v2c, default=10000) / 2, max(table['sum_cov'] * 0.95, default=10000)), y=y_cut,
                                 va='bottom', ha='center', s=str(value))

                    otab.write('\t'.join(map(str, vals)) + '\n')
                    plt.savefig(out_folder + cor_folder + '_'.join([cell_line, model.replace('.', ','), 'cor']))

                if 'cor-seg' in plots:
                    fig = plt.figure(figsize=(10, 8))
                    sns.scatterplot(data=table,
                                    x='number_of_segments_' + model,
                                    y='cor_by_snp_' + model,
                                    color='green')
                    plt.vlines(x=np.mean(table['total_regions']), ymin=-1, ymax=1)
                    plt.hlines(y=0.3, xmin=0, xmax=max(table['number_of_segments_' + model]))
                    value = len(table[table['cor_by_snp_' + model] >= 0.3].index)
                    plt.text(x=max(table['number_of_segments_' + model]) / 2, y=0.3, va='bottom', s=str(value))
                    plt.hlines(y=0, xmin=0, xmax=max(table['number_of_segments_' + model]))
                    value = len(table[table['cor_by_snp_' + model] >= 0].index)
                    plt.text(x=max(table['number_of_segments_' + model]) / 2, y=0, va='bottom', s=str(value))
                    plt.title('\n'.join((cell_line, model, 'Correlation vs # of segments', 'NEW N')))
                    plt.grid(True)
                    plt.savefig(out_folder + cor_seg_folder + '_'.join([cell_line, model.replace('.', ','), 'cor-seg']))

                if 'cosm_heatmap' in plots:
                    heat_table = None
                    f1, f2 = regs[heatmap_reg]
                    if f2 is None: f2 = max(table['total_snps'])
                    f1, f2 = map(float, (f1, f2))
                    for file_name in os.listdir(heatmapdata_path + model + '_tables/'):
                        cell_line_name = file_name[:file_name.rfind('_')]
                        if cell_line_name == cell_line:
                            try:
                                t = pd.read_csv(heatmapdata_path + model + '_tables/' + file_name, sep='\t')
                            except Exception as e:
                                print(e.args[0])
                                continue
                            t.columns = ['chr', 'pos', 'value', 'segment_value']
                            t = t[t['value'] != 0.0][['value', 'segment_value']]
                            print(file_name)
                            r, _ = t.shape
                            if r < f1 or r >= f2:
                                print('skipping, r={}'.format(r))
                                continue
                            if heat_table is None:
                                heat_table = t
                            else:
                                heat_table = heat_table.append(t)

                    heat_table = heat_table.rename(columns={'segment_value': 'COSMIC', 'value': 'SEGMENTATION'})
                    heat_matrix = pd.crosstab(heat_table['SEGMENTATION'], heat_table['COSMIC'])
                    # heat_table = None
                    heat_matrix.sort_index(ascending=False, inplace=True)
                    order = heat_matrix.index.tolist()
                    rename_dict = {1.0: '1', 1.3333333333333337: '4/3', 1.5: '1.5', 2.0: '2', 2.5: '5/2', 3.0: '3',
                                   3.666666666666667: '11/3', 4.0: '4', 5.0: '5', 6.0: '6', 7.0: '7', 8.0: '8',
                                   9.0: '9', 10.0: '10', 11.0: '11'}
                    heat_matrix = heat_matrix.rename(index=rename_dict, columns=rename_dict)
                    cosmic_dist = heat_matrix.sum(axis=0)
                    seg_dist = heat_matrix.sum(axis=1)
                    heat_freq = heat_matrix.div(cosmic_dist, axis=1)
                    heat_freq = heat_freq * 100

                    fs = (max(1.3 * len(heat_matrix.columns), 8), max(1.2 * len(heat_matrix.index), 6))
                    print(fs)
                    fig = plt.figure(tight_layout=True, figsize=fs)
                    gs = fig.add_gridspec(5, 10)
                    ax1 = fig.add_subplot(gs[0, 1:-2])
                    ax2 = fig.add_subplot(gs[1:, 1:-2])
                    ax3 = fig.add_subplot(gs[1:, -2:])
                    ax_c = fig.add_subplot(gs[1:, 0])

                    sns.heatmap(heat_freq, annot=heat_matrix, fmt='d', vmin=0, vmax=100,
                                cmap=sns.cubehelix_palette(200, start=.5, rot=-.75, reverse=True),
                                cbar=True, ax=ax2, cbar_ax=ax_c)
                    ax_c.yaxis.set_label_position('left')
                    ax_c.yaxis.set_ticks_position('left')
                    sns.countplot(ax=ax1, x='COSMIC', data=heat_table, color='grey')
                    ax1.set(xlabel='', ylabel='', xticklabels=[], xticks=[], yticklabels=[], yticks=[])
                    sns.countplot(ax=ax3, y='SEGMENTATION', data=heat_table, color='grey', order=order)
                    ax3.set(xlabel='', ylabel='', xticklabels=[], xticks=[], yticklabels=[], yticks=[])
                    plt.suptitle('\n'.join([short_name[cell_line], '\n'.join(model.split('-')), regnames[heatmap_reg]]),
                                 ha='left', x=9 / 11,
                                 y=0.9)
                    plt.savefig(out_folder + heat_folder + regnames[heatmap_reg] + '/' + '_'.join(
                        [cell_line, model.replace('.', ','), 'heatmap']))
