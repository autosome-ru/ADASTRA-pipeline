import os.path
import pandas as pd


from scripts.HELPERS.paths import create_path_from_master_list_df, get_result_stats_path
from scripts.HELPERS.paths_for_components import master_list_path
from scripts.HELPERS.helpers import dtype_dict, states


def collect_stats(key_names=None, BAD=None, suffix=''):
    ml_df = pd.read_table(master_list_path, dtype=dtype_dict)
    out_t = None

    for index, row in ml_df.iterrows():
        if key_names is not None:
            if row['CELLS'] not in key_names:  # <------
                continue
        bad_table_path = create_path_from_master_list_df(row, 'BAD')
        if not os.path.isfile(bad_table_path):
            continue
        df = pd.read_table(bad_table_path)
        if df.empty:
            continue
        if BAD is not None:
            sum_df = df[df['BAD'] == BAD][['ref_read_counts', 'alt_read_counts']]  # <------
        else:
            sum_df = df[['ref_read_counts', 'alt_read_counts']]

        if out_t is None:
            out_t = pd.DataFrame()
            out_t['alt_counts'] = sum_df['alt_read_counts']
            out_t['ref_counts'] = sum_df['ref_read_counts']
            out_t = out_t.groupby(['alt_counts', 'ref_counts']).size().reset_index(name='counts')
            out_t.fillna(0, inplace=True)
        else:
            tmp_df = pd.DataFrame()
            tmp_df['alt_counts'] = sum_df['alt_read_counts']
            tmp_df['ref_counts'] = sum_df['ref_read_counts']
            tmp_df = tmp_df.groupby(['alt_counts', 'ref_counts']).size().reset_index(name='counts')
            tmp_df.fillna(0, inplace=True)
            out_t = out_t.append(tmp_df).groupby(['alt_counts', 'ref_counts'], as_index=False).sum()
    if out_t is None:
        print('No stats collected')
        return
    out_t.to_csv(os.path.join(get_result_stats_path(),
                              'fixed_alt_bias_statistics_BAD={:.1f}{}.tsv'.format(BAD if BAD else 0, suffix)),
                 sep="\t", index=False)


def main():
    for BAD in states:
        collect_stats(None, BAD=BAD)


# if __name__ == "__main__":
#     esc = {'H1 (embryonic stem cells)',
#              'BGO3 (embryonic stem cells)',
#              'HUES64 (embryonic stem cells)',
#              'CyT49 (embryonic stem cells)',
#              'H9 (embryonic stem cells)',
#              'VAL-3 (embryonic stem cells)',
#              'WA09 (embryonic stem cells)',
#              'embryonic stem cells',
#              'human embryonic stem cells, H1 (WA01)'}
#     hct116 = {"HCT-116 (colon carcinoma)"}
#     k562 = {'K562 (myelogenous leukemia)'}
#     for BAD in [None, 1, 2, 3, 4, 5, 6, 4/3, 5/2, 3/2]:
#         collectFixedAltStatistics(BAD=BAD, key_name=None, suffix='')
#         collectFixedAltStatistics(BAD=BAD, key_name=esc, suffix='_esc')
#         collectFixedAltStatistics(BAD=BAD, key_name=hct116, suffix='_hct116')
#         collectFixedAltStatistics(BAD=BAD, key_name=k562, suffix='_k562')
#     main()
