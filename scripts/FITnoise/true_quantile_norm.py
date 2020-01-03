import os
import sys
import numpy as np
import pandas as pd


def get_normalized_allele_reads(df):
    # print(np.cumsum(df['ref']), np.cumsum(df['alt']))
    prev = 0
    cur_mean_num = 0
    cur_mean_denom = 0
    cur_ar_index = 0
    # cur_ar = df['allele_reads'][cur_ar_index]
    cur_new_ar_index = 0
    new_allele_reads = np.zeros(len(df.index), dtype=np.float_)
    iterator = sorted(sorted(list(zip(np.cumsum(df['alt']), ['alt'] * len(df.index))) + \
                             list(zip(np.cumsum(df['ref']), ['ref'] * len(df.index))), key=lambda x: x[1],
                             reverse=True),
                      key=lambda x: x[0])
    for cs, marker in iterator:
        cur_mean_num += (cs - prev) * df['allele_reads'][cur_ar_index]
        cur_mean_denom += cs - prev

        if marker == 'ref':
            if cur_mean_denom != 0:
                new_allele_reads[cur_new_ar_index] = cur_mean_num / cur_mean_denom
            else:
                new_allele_reads[cur_new_ar_index] = df['allele_reads'][cur_ar_index]
            cur_new_ar_index += 1
            cur_mean_num = 0
            cur_mean_denom = 0

        elif marker == 'alt':  # explicit better than implicit
            cur_ar_index += 1

        prev = cs
    return new_allele_reads


if __name__ == '__main__':
    for BAD in [1, 2, 3, 4, 5, 6, 4 / 3, 5 / 2, 3 / 2]:
        stats = pd.read_table(os.path.expanduser('~/bias_statistics_BAD={:.1f}.tsv'.format(BAD)))
        stats['new_allele_reads'] = get_normalized_allele_reads(stats)
        print(stats)
        stats[['allele_reads', 'new_allele_reads']].to_csv(
            os.path.expanduser('~/ref_counts_scaling_BAD={:.1f}.tsv'.format(BAD)), index=False, sep='\t')
