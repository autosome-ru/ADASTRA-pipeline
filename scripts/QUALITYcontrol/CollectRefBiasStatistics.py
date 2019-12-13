import sys
import os.path
import json
import pandas as pd

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import cl_dict_path, parameters_path


def collectRefAltStatistics(key_name=None, BAD=None):
    with open(cl_dict_path, "r") as read_file:
        cell_lines_dict = json.loads(read_file.readline())
    out_t = None

    for key in cell_lines_dict:
        if key_name is not None:
            if key not in key_name:  # <------
                continue
        for align_path in cell_lines_dict[key]:
            if not os.path.isfile(align_path):
                continue
            print(align_path)
            df = pd.read_table(align_path)
            if df.empty:
                continue
            if BAD is not None:
                sum_df = df[df['BAD'] == BAD][['ref_read_counts', 'alt_read_counts']]  # <------
            else:
                sum_df = df[['ref_read_counts', 'alt_read_counts']]

            if out_t is None:
                out_t = pd.DataFrame()
                out_t['ref'] = sum_df['ref_read_counts'].value_counts()
                out_t['alt'] = sum_df['alt_read_counts'].value_counts()
                out_t['allele_reads'] = out_t.index
                out_t = out_t.reset_index(drop=True)
                out_t.fillna(0, inplace=True)
            else:
                tmp_df = pd.DataFrame()
                tmp_df['ref'] = sum_df['ref_read_counts'].value_counts()
                tmp_df['alt'] = sum_df['alt_read_counts'].value_counts()
                tmp_df['allele_reads'] = tmp_df.index
                tmp_df = tmp_df.reset_index(drop=True)
                tmp_df.fillna(0, inplace=True)
                out_t = out_t.append(tmp_df).groupby('allele_reads', as_index=False).sum()
    if out_t is None:
        return
    with open(parameters_path + 'bias_statistics.tsv', 'w') as out:
        out_t.to_csv(out, sep="\t")


def collectCoverStatistics(key_name=None, BAD=None):
    with open(cl_dict_path, "r") as read_file:
        cell_lines_dict = json.loads(read_file.readline())
    out_t = None

    for key in cell_lines_dict:
        if key_name is not None:
            if key not in key_name:  # <------
                continue
        for align_path in cell_lines_dict[key]:
            if not os.path.isfile(align_path):
                continue
            print(align_path)
            df = pd.read_table(align_path)
            if df.empty:
                continue
            if BAD is not None:
                sum_df = df[df['BAD'] == BAD][['ref_read_counts', 'alt_read_counts']]  # <------
            else:
                sum_df = df[['ref_read_counts', 'alt_read_counts']]

            if out_t is None:
                out_t = pd.DataFrame()
                out_t['cover'] = (sum_df['ref_read_counts'] + sum_df['alt_read_counts'])
                out_t['ref_counts'] = sum_df['ref_read_counts']
                out_t = out_t.groupby(['cover', 'ref_counts']).size().reset_index(name='counts')
                out_t.fillna(0, inplace=True)
            else:
                tmp_df = pd.DataFrame()
                tmp_df['cover'] = (sum_df['ref_read_counts'] + sum_df['alt_read_counts'])
                tmp_df['ref_counts'] = sum_df['ref_read_counts']
                tmp_df = tmp_df.groupby(['cover', 'ref_counts']).size().reset_index(name='counts')
                tmp_df.fillna(0, inplace=True)
                out_t = out_t.append(tmp_df).groupby(['cover', 'ref_counts'], as_index=False).sum()
    if out_t is None:
        return
    with open(parameters_path + 'cover_bias_statistics.tsv', 'w') as out:
        out_t.to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    embryonic = ['H1_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'BGO3_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'H9_embryonic_stem_cells',
                 'WA09_embryonic_stem_cells',
                 'human_embryonic_stem_cells,_H1_WA01',
                 'H1_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'HUES64_embryonic_stem_cells',
                 'H9_embryonic_stem_cells',
                 'H9_embryonic_stem_cells',
                 'H9_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'embryonic_stem_cells',
                 'CyT49_embryonic_stem_cells',
                 'H9_embryonic_stem_cells',
                 'H9_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'VAL-3_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'embryonic_stem_cells',
                 'embryonic_stem_cells',
                 'H9_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'BGO3_embryonic_stem_cells',
                 'H9_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'H9_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'embryonic_stem_cells',
                 'H9_embryonic_stem_cells',
                 'embryonic_stem_cell-derived_pancreatic_cells',
                 'H9_embryonic_stem_cells',
                 'embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'HUES64_embryonic_stem_cells',
                 'H9_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'H1_embryonic_stem_cells',
                 'embryonic_stem_cell-derived_neural_progenitors',
                 'embryonic_stem_cell-_derived_mid_radial_glial_neural_progenitors',
                 'H9_embryonic_stem_cells']
    all_diploids = ['226LDM (normal breast luminal cells)',
                    'AB32 (normal breast cells)',
                    'adrenal gland',
                    'AG04449 (fibroblasts)', 'AG04450 (fibroblasts)',
                    'AG09309 (fibroblasts)', 'AG09319 (fibroblasts)', 'AG10803 (fibroblasts)',
                    'AM-iPS-6 (induced pluripotent cells)', 'aortic adventitial fibroblasts', 'ascending aorta',
                    'astrocytes', 'B-cells', 'BEAS-2B (bronchial epithelium)', 'bipolar neuron', 'BJ (fibroblasts)',
                    'blood monocytes', 'breast epithelial cells', 'bronchial epithelial cells', 'CD14+ monocytes',
                    'CD34+ hematopoietic stem cells', 'CD34+ hematopoietic stem cells-derived proerythroblasts',
                    'CD34+ hematopoietic stem progenitor cells', 'CD34+ stem cells-derived erythroblasts',
                    'CD36+ erythroid cells', 'CD4+ T-cells', 'CD8+ T-cells', 'cord blood-derived mononuclear cells',
                    'coronary artery', 'coronary artery smooth muscle', 'cranial neural crest cells',
                    'CyT49 (embryonic stem cells)', 'CyT49-derived endodermal cells', 'dermal fibroblasts',
                    'E14 retina', 'E16 retina', 'E18 retina', 'E23 retina',
                    'Edom-iPS-2 (induced pluripotent stem cells)', 'embryonic kidney', 'embryonic kidney cortex',
                    'embryonic stem cells', 'endometrial stromal cells', 'epidermal keratinocytes',
                    'epidermal stem cells', 'epithelial cell of prostate', 'epithelial cell of proximal tubule',
                    'erythroblasts', 'erythroid cells', 'erythroid progenitors', 'esophageal epithelial cells',
                    'esophagus muscularis mucosa', 'esophagus squamous epithelium', 'FB0167P (progeria fibroblasts)',
                    'FB8470 (fibroblasts)', 'fetal CD34+ hematopoietic stem progenitor cells', 'fetal lung fibroblasts',
                    'fetal osteoblasts', 'fetal proerythroblasts', 'fetal prostate fibroblasts',
                    'fibroblast of dermis', 'fibroblast of mammary gland', 'fibroblast of pulmonary artery',
                    'Flp-In-293 (embryonic kidney)', 'foreskin', 'foreskin keratinocyte', 'gastrocnemius medialis',
                    'gastroesophageal sphincter', 'germinal center B-cells', 'GM06990 (female B-cells)',
                    'GM10248 (male B-cells)', 'GM10266 (male B-cells)', 'GM10855 (female B-cells)',
                    'GM10861 (female B-cells)', 'GM12864 (male B-cells)', 'GM12865 (female B-cells)',
                    'GM12866 (male B-cells)', 'GM12867 (female B-cells)', 'GM12868 (female B-cells)',
                    'GM12869 (female B-cells)', 'GM12870 (male B-cells)', 'GM12871 (male B-cells)',
                    'GM12872 (male B-cells)', 'GM12873 (female B-cells)', 'GM12874 (male B-cells)',
                    'GM12875 (female B-cells)', 'GM12878 (female B-cells)', 'GM12891 (male B-cells)',
                    'GM12892 (female B-cells)', 'GM13976 (female B-cells)', 'GM13977 (female B-cells)',
                    'GM15510 (B-cells)', 'GM18505 (female B-cells)', 'GM19099 (female B-cells)',
                    'GM19238 (female B-cells)', 'GM19239 (male B-cells)', 'GM19240 (B-cells)',
                    'GM20000 (B-cells)', 'GM23248', 'GM23338', 'H1 (embryonic stem cells)',
                    'H1-derived mesenchymal stem cells', 'H1-derived mesendodermal cells',
                    'H1-derived neuronal progenitors', 'H1-derived trophectodermal cells',
                    'H9 (embryonic stem cells)', 'H9-derived neural progenitors',
                    'H9-derived neuroectoderm cells', 'HAEC (human aortic endothelial cells)', 'hASC (preadipocytes)',
                    'HCPEpiC (choroid plexus epithelial cells)', 'heart left ventricle', 'HEK293 (embryonic kidney)',
                    'HEK293T (embryonic kidney)', 'hematopoietic stem and progenitor cells',
                    'hematopoietic stem cells and progenitors', 'HFF-Myc (fetal fibroblasts)',
                    'hMADS-3 (adipose-derived stem cells)', 'HUES64 (embryonic stem cells)',
                    'human astrocytes-cerebellar', 'human astrocytes-spinal cord',
                    'human brain microvascular endothelial cells', 'human cardiac fibroblasts-adult atrial',
                    'Human Umbilical Cord Blood-Derived Erythroid Progenitor Cells-2 (HUDEP-2)',
                    'human villous mesenchymal fibroblast cells', 'HUVEC (umbilical vein endothelial cells)',
                    'IB4 (lymphocytes)', 'IMR90 (lung fibroblasts)', 'induced pluripotent stem cells',
                    'IPSC-derived bipolar neuron', 'keratinocytes', 'kidney', 'kidney epithelial cell',
                    'lower leg skin', 'lung', 'lymphoblastoid cells', 'macrophages', 'mammary epithelial cells',
                    'MCF10A (breast epithelial cells)', 'mesenchymal stem cells', 'monocyte-derived dendritic cells',
                    'monocyte-derived macrophages', 'MRC-5 (fetal lung fibroblasts)',
                    'MRC-iPS-25 (induced pluripotent stem cells)', 'myoblasts', 'myofibroblasts',
                    'neonatal keratinocytes', 'neural cell', 'neural progenitors', 'neural stem cells', 'neutrophils',
                    'normal epidermal keratinocytes', 'normal lung fibroblasts', 'omental fat pad', 'osteoblasts',
                    'ovary', 'pancreas', 'patient-derived pleaural effusion', 'periferal blood macrophages',
                    'periferal blood stem cells', 'peripheral blood mononuclear cells', "Peyer's patch",
                    'PrEC (prostate cells)', 'primary endothelial colony-forming cells', 'primary erythroblasts',
                    'primary foreskin fibroblasts', 'proerythroblasts', 'prostate gland', 'renal cortex',
                    'renal tubular cells', 'retinal pigment epithelial cells', 'right atrium auricular region',
                    'RWPE1 (prostate epithelial cells)', 'sigmoid colon', 'skeletal muscles and myoblasts',
                    'smooth muscle cell', 'spleen', 'stem and progenitor cells', 'stomach',
                    'subcutaneous adipose tissue', 'subcutaneous white adipose tissue', 'suprapubic skin',
                    'T-cells', 'testis', 'Th1-cells', 'Th2-cells', 'thoracic aorta', 'thymus',
                    'thymus-derived CD34- cells', 'thymus-derived CD34+ cells', 'thyroid gland',
                    'tibial artery', 'tibial nerve', 'tracheobronchial epithelial cells', 'transverse colon',
                    'T-REx-Jurkat (T-cells)', 'upper lobe of left lung', 'UtE-iPS-4 (induced pluripotent stem cells)',
                    'UtE-iPS-6 (induced pluripotent stem cells)', 'UtE-iPS-7 (induced pluripotent stem cells)',
                    'uterus', 'vagina', 'WI-38 (lung fibroblasts)']
    collectCoverStatistics(key_name=all_diploids, BAD=4/3)
