import os
import re
import sys
import json
from scipy import stats as st
import numpy as np
import pandas as pd

from scripts.HELPERS.paths_for_components import badmaps_path, badmaps_dict_path, correlation_path
from scripts.HELPERS.helpers import Intersection, pack, UnpackBadSegments, get_states


def unpack_snps(line):
    line = line.strip().split('\t')
    return [line[0], int(line[1]), int(line[5]), int(line[6]), line[7]]


def get_p_value(n, p, x):
    dist = st.binom(n=n, p=p)
    cdf = dist.cdf
    return (cdf(x) - cdf(4) + cdf(n-5) - cdf(n-x-1)) / (cdf(n-5) - cdf(4)) if x < n/2 else 1


def main(file_name):
    if not os.path.isdir(correlation_path):
        try:
            os.mkdir(correlation_path)
        except:
            pass

    with open(badmaps_dict_path, 'r') as file:
        aligns_by_cell_type = json.loads(file.readline().strip())

    modes = []
    for dir_name in sorted(os.listdir(badmaps_path)):
        if os.path.isdir(os.path.join(badmaps_path, dir_name)) and dir_name != 'merged_vcfs':
            modes.append(dir_name)

    try:
        assert os.path.isfile(os.path.join(badmaps_path, 'merged_vcfs', file_name))
    except AssertionError:
        print(os.path.join(badmaps_path, 'merged_vcfs', file_name), file_name)
        exit(1)

    name = file_name.split('@')[0]
    lab = file_name.split('@')[1][:-4]  # .tsv

    try:
        aligns = aligns_by_cell_type[file_name[:-4]]  # .tsv
        al_list = [os.path.basename(align) for align in aligns if os.path.isfile(align)]
        datasetsn = len(al_list)
    except KeyError:
        datasetsn = 'nan'
        al_list = []
        print(file_name)

    table_path = os.path.join(badmaps_path, 'merged_vcfs', file_name)
    for mode in modes:
        if re.match(r'^CAIC@.+@.+$', mode) is not None:
            states = get_states(mode.split('@')[1])
        else:
            states = get_states('')
        if not os.path.isdir(os.path.join(correlation_path, mode + '_tables')):
            try:
                os.mkdir(os.path.join(correlation_path, mode + '_tables'))
            except:
                pass
        badmaps_file_path = os.path.join(badmaps_path, mode, name + '@' + lab + '.badmap.tsv')
        out_path = os.path.join(correlation_path, mode + '_tables', name + '_' + lab.replace('_', '-') + '.tsv')
        print(out_path)

        u = UnpackBadSegments(0)

        with open(table_path, 'r') as table, open(badmaps_file_path, 'r') as BADmap_file, open(out_path, 'w') as out:
            out.write('#' + str(datasetsn) + '@' + lab + '@' + ','.join(al_list) + '\n')
            for chr, pos, ref, alt, filename, in_intersection, segment_BAD, segment_id, Qual, segn, sumcov \
                    in Intersection(table, BADmap_file,
                                    unpack_segments_function=lambda x: u.unpackBADSegments(x, states), unpack_snp_function=unpack_snps,
                                    write_intersect=True, write_segment_args=True):
                if not in_intersection:
                    continue
                p_value = get_p_value(ref + alt, 1 / (segment_BAD + 1), min(ref, alt))
                out.write(pack([chr, pos, ref, alt, segment_BAD] + [Qual[x] for x in Qual] + [segn, sumcov] + [filename, segment_id, p_value]))


if __name__ == '__main__':
    redo = ['HEK293__embryonic_kidney_@_labs_michael-snyder___biosamples_ENCBS957THJ_.tsv', '22RV1__prostate_carcinoma_@GSE94013.tsv', 'LNCaP__prostate_carcinoma_@GSE72467.tsv', 'E18_retina@GSE86981.tsv', 'primary_monocytes@GSE136216.tsv', 'hES-PRE@GSE114305.tsv', 'LS180__colon_cancer_@GSE73319.tsv', 'K562__myelogenous_leukemia_@_labs_michael-snyder___biosamples_ENCBS041GBQ_.tsv', 'MCF-7_C4-12__invasive_breast_ductal_carcinoma_@GSE48096.tsv', 'stomach@_labs_michael-snyder___biosamples_ENCBS424ANS_.tsv', 'HMEC__human_mammary_epithelial_cells_@GSE62425.tsv', 'GM06990__female_B-cells_@_labs_john-stamatoyannopoulos___biosamples_ENCBS185XLT_.tsv', 'T3M-1_Cl-10__oral_cavity_squamous_cell_carcinoma_@GSE53357.tsv', 'prostate_cancer-associated_fibroblasts@GSE90772.tsv', 'gastroesophageal_sphincter@_labs_richard-myers___biosamples_ENCBS057FSI_.tsv', 'A549__lung_carcinoma_@_labs_tim-reddy___biosamples_ENCBS293HYY_.tsv', 'HeLa_S3__cervical_adenocarcinoma_@GSE25416.tsv', 'K562__myelogenous_leukemia_@_labs_michael-snyder___biosamples_ENCBS398HUM_.tsv', 'HepG2__hepatoblastoma_@_labs_xiang-dong-fu___biosamples_ENCBS604BAZ_.tsv', 'K562__myelogenous_leukemia_@GSE24777.tsv', 'MCF7__Invasive_ductal_breast_carcinoma_@GSE46055.tsv', 'HEK293__embryonic_kidney_@GSE97540.tsv', 'foreskin_keratinocyte@_labs_michael-snyder___biosamples_ENCBS079LRU_.tsv', 'epidermal_stem_cells@GSE65838.tsv', 'BJAB__Burkitt_lymphoma_@GSE86154.tsv', 'subcutaneous_adipose_tissue@_labs_richard-myers___biosamples_ENCBS647IID_.tsv', 'LNCaP__prostate_carcinoma_@GSE52201.tsv', 'HepG2__hepatoblastoma_@_labs_peggy-farnham___biosamples_ENCBS950MSQ_.tsv', 'LNCaP__prostate_carcinoma_@GSE48308.tsv', 'HepG2__hepatoblastoma_@_labs_richard-myers___biosamples_ENCBS543FMH_.tsv', 'HEK293T__embryonic_kidney_@GSE62616.tsv', 'MCF7__Invasive_ductal_breast_carcinoma_@GSE44737.tsv', 'GM12878__female_B-cells_lymphoblastoid_cell_line_@_labs_michael-snyder___biosamples_ENCBS389LEA_.tsv', 'AGS__gastric_adenocarcinoma_@GSE51936.tsv', 'MCF7__Invasive_ductal_breast_carcinoma_@GSE74141.tsv', 'A549__lung_carcinoma_@_labs_tim-reddy___biosamples_ENCBS393PAL_.tsv', 'Dombi23__primary_keratinocytes_@GSE123711.tsv', 'HepG2__hepatoblastoma_@_labs_richard-myers___biosamples_ENCBS980ARE_.tsv', 'HEL__Erythro_Leukemia_@GSE112265.tsv', 'HEK293__embryonic_kidney_@GSE68714.tsv', 'HepG2__hepatoblastoma_@_labs_michael-snyder___biosamples_ENCBS157PJQ_.tsv', 'MDA-MB-453__breast_adenocarcinoma_@GSE45201.tsv', 'HCT-116__colon_carcinoma_@GSE57398.tsv']
    if sys.argv[1] not in redo:
        sys.exit(0)
    main(sys.argv[1])
