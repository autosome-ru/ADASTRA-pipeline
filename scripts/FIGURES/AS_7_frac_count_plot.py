import sys
import numpy as np
import os
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


def get_color(row):
    if abs(row['log_fc']) < fc_tr:
        return 'N'
    # if row[field + '_ref'] < fdr_tr and row[field + '_alt'] < fdr_tr:
    #     return 'purple'
    if row['log_fc'] * row['log_pv'] > 0:
        return 'C'
    else:
        return 'D'


if __name__ == '__main__':
    FDR = 1
    ES = 1
    TF = 1

    sns.set(font_scale=1.4, style="ticks", font="lato", palette=('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2',
                                                                 '#D55E00', '#CC79A7'))
    sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})
    plt.rcParams['font.weight'] = "medium"
    plt.rcParams['axes.labelweight'] = 'medium'
    plt.rcParams['figure.titleweight'] = 'medium'
    plt.rcParams['axes.titleweight'] = 'medium'
    plt.rcParams['figure.figsize'] = 7, 8
    plt.rcParams["legend.framealpha"] = 1
    plt.rcParams['axes.xmargin'] = 0
    plt.rcParams['axes.ymargin'] = 0
    plt.rcParams["legend.framealpha"] = 1

    field = 'fdrp_bh'

    perf_tr = 0.0005
    fc_tr = 1

    title = 'BillCipher'

    # FDR
    if FDR:
        pt = pd.read_table(os.path.expanduser("C:\\Users\\boyts\OneDrive\Desktop/BillCipherStats/scatter/{}.tsv".format('all_tfs')))

        pt = pt[(pt['motif_log_pref'] >= -np.log10(perf_tr)) | (pt['motif_log_palt'] >= -np.log10(perf_tr))]

        pt = pt[~(pt[field + '_alt'].isnull() | pt[field + '_ref'].isnull())]
        pt = pt[~(pt['motif_fc'].isnull())]
        pt['log_pv'] = (np.log10(
            pt[[field + '_ref', field + '_alt']]).min(axis=1)) \
                       * np.sign(pt[field + '_alt'] - pt[field + '_ref'])
        pt['log_fc'] = pt['motif_fc']

        pt['col'] = pt.apply(lambda x: get_color(x), axis=1)

        grid = [0.1 * x for x in range(301)]
        fraction = []
        total_number_of_concordant_and_discordant_snps_at_treshold = []

        for cutoff in grid:
            # print(cutoff)
            pt = pt[(pt['log_pv'] >= cutoff) | (pt['log_pv'] <= -cutoff)]
            blue = len(pt[pt['col'] == 'C'].index)
            red = len(pt[pt['col'] == 'D'].index)
            grey = len(pt[pt['col'] == 'N'].index)
            fraction.append(blue/(blue+red))
            total_number_of_concordant_and_discordant_snps_at_treshold.append(blue + red)

        # ct = 5000
        # c_fdr = 0
        # index = 0
        # for ind, (cut, tot) in enumerate(zip(grid, total_number_of_concordant_and_discordant_snps_at_treshold)):
        #     if tot >= ct:
        #         c_fdr = cut
        #         index = ind
        # print(title)
        # print(c_fdr, fraction[index])

        fig, ax1 = plt.subplots()
        plt.tight_layout(pad=2.5)

        ax1.plot(total_number_of_concordant_and_discordant_snps_at_treshold, fraction, zorder=2)
        ax1.grid(True)
        ax1.set_ylim(0.5, 1)
        ax1.scatter(x=total_number_of_concordant_and_discordant_snps_at_treshold[13],
                    y=fraction[13],
                    color='grey')
        ax1.scatter(x=total_number_of_concordant_and_discordant_snps_at_treshold[20],
                    y=fraction[20],
                    color='grey')
        ax1.set_ylabel('Fraction of concordant SNPs')
        ax1.set_yticks([0.6, 0.7, 0.8, 0.9, 1])
        ax1.set_xscale('log')
        ax1.set_xlabel('Number of concordant\nand discordant SNPs')
        ax1.set_xlim(10, 100000)

        plt.suptitle(title)

        plt.savefig(os.path.expanduser("C:\\Users\\boyts\OneDrive\Desktop/BillCipherStats/AS_7/fraction_count.png"), dpi=300)
        plt.show()
        plt.close(fig)

    # ES
    if ES:
        pt = pd.read_table(os.path.expanduser("C:\\Users\\boyts\OneDrive\Desktop/BillCipherStats/scatter/{}.tsv".format('all_tfs')))

        pt = pt[(pt['motif_log_pref'] >= -np.log10(perf_tr)) | (pt['motif_log_palt'] >= -np.log10(perf_tr))]

        pt = pt[~(pt[field + '_alt'].isnull() | pt[field + '_ref'].isnull())]
        pt = pt[~(pt['motif_fc'].isnull())]
        pt['log_pv'] = (np.log10(
            pt[[field + '_ref', field + '_alt']]).min(axis=1)) \
                       * np.sign(pt[field + '_alt'] - pt[field + '_ref'])
        pt['log_fc'] = pt['motif_fc']
        pt['log2_es'] = pt[['es_mean_ref', 'es_mean_alt']].max(axis=1) / np.log(2)

        pt['col'] = pt.apply(lambda x: get_color(x), axis=1)

        grid = [0.1 * x for x in range(35)]
        fraction = []
        total_number_of_concordant_and_discordant_snps_at_treshold = []

        for cutoff in grid:
            # print(cutoff)
            pt = pt[(pt['log2_es'] >= cutoff)]
            blue = len(pt[pt['col'] == 'C'].index)
            red = len(pt[pt['col'] == 'D'].index)
            grey = len(pt[pt['col'] == 'N'].index)
            fraction.append(blue/(blue+red))
            total_number_of_concordant_and_discordant_snps_at_treshold.append(blue + red)

        fig, (ax1, ax2) = plt.subplots(2, 1)
        plt.tight_layout(pad=2.5)

        ax1.plot(grid, fraction)
        ax1.grid(True)
        ax1.set_ylim(0.5, 1)
        ax1.set_ylabel('Fraction of concordant SNPs')
        ax1.set_yticks([0.6, 0.7, 0.8, 0.9, 1])
        ax2.plot(grid, total_number_of_concordant_and_discordant_snps_at_treshold)
        ax2.set_yscale('log')
        ax2.grid(True)
        ax2.set_xlabel('log₂ effect size')
        ax2.set_ylabel('Number of concordant\nand discordant SNPs')

        plt.suptitle(title)

        plt.savefig(os.path.expanduser("C:\\Users\\boyts\OneDrive\Desktop/BillCipherStats/AS_7/fraction_count_es.png"), dpi=300)
        plt.show()
        plt.close(fig)

    #TF
    if TF:
        for tf_name in 'CTCF_HUMAN', 'ANDR_HUMAN':
            pt = pd.read_table(os.path.expanduser("C:\\users\\boyts\\OneDrive\Desktop\BillCipherStats/scatter/{}.tsv").format(tf_name))

            pt = pt[(pt['motif_log_pref'] >= -np.log10(perf_tr)) | (pt['motif_log_palt'] >= -np.log10(perf_tr))]

            pt = pt[~(pt[field + '_alt'].isnull() | pt[field + '_ref'].isnull())]
            pt = pt[~(pt['motif_fc'].isnull())]
            pt['log_pv'] = (np.log10(
                pt[[field + '_ref', field + '_alt']]).min(axis=1)) \
                           * np.sign(pt[field + '_alt'] - pt[field + '_ref'])
            pt['log_fc'] = pt['motif_fc']

            pt['col'] = pt.apply(lambda x: get_color(x), axis=1)

            grid = [0.01 * x for x in range(2001)]
            fraction = []
            total_number_of_concordant_and_discordant_snps_at_treshold = []

            for cutoff in grid:
                # print(cutoff)
                pt = pt[(pt['log_pv'] >= cutoff) | (pt['log_pv'] <= -cutoff)]
                blue = len(pt[pt['col'] == 'C'].index)
                red = len(pt[pt['col'] == 'D'].index)
                grey = len(pt[pt['col'] == 'N'].index)
                fraction.append(blue / (blue + red))
                total_number_of_concordant_and_discordant_snps_at_treshold.append(blue + red)

            if tf_name == 'CTCF_HUMAN':
                mat_fdr = 0
                for cut, tot in zip(grid, total_number_of_concordant_and_discordant_snps_at_treshold):
                    if tot >= 505 + 47:
                        mat_fdr = cut
                print(mat_fdr)

            fig, ax1 = plt.subplots()
            plt.tight_layout(pad=2.5)

            ax1.plot(total_number_of_concordant_and_discordant_snps_at_treshold, fraction, zorder=2)
            ax1.grid(True)
            ax1.set_ylim(0.5, 1)
            ax1.scatter(x=total_number_of_concordant_and_discordant_snps_at_treshold[130],
                        y=fraction[130],
                        color='grey')
            ax1.scatter(x=total_number_of_concordant_and_discordant_snps_at_treshold[200],
                        y=fraction[200],
                        color='grey')
            ax1.set_ylabel('Fraction of concordant SNPs')
            ax1.set_yticks([0.6, 0.7, 0.8, 0.9, 1])
            ax1.set_xscale('log')
            ax1.set_xlabel('Number of concordant\nand discordant SNPs')
            ax1.set_xlim(10, 10000)
            if tf_name == 'CTCF_HUMAN':
                ax1.scatter(505 + 47, [505/(505+47)], color='C2', zorder=3)
                ax1.text(x=505 + 47, y=505/(505+47), color='C2', va='bottom', ha='left', s='Shi et al.')

            plt.suptitle(title)

            plt.savefig(os.path.expanduser("C:\\Users\\boyts\OneDrive\Desktop/BillCipherStats/AS_7/fraction_count_{}.png".format(tf_name)), dpi=300)
            plt.show()
            plt.close(fig)