import glob
import os
from .annotate_with_phenotypes import main as annotate_main, phenotype_db_names
from scripts.HELPERS.paths import get_result_dir_path, get_release_stats_path
from scripts.HELPERS.helpers import check_and_create_dir
import shutil


def Get_ASB(mode, release_path, pv='fdrp_bh_', fdr=0.05, peaks='>=', nagg=''):
    snps = {}

    files = glob.glob('{}/*.tsv'.format(release_path))

    for tsv in sorted(files):

        count = 0

        tf = {'TF': tsv[tsv.rfind('/') + 1:tsv.rfind('_HUMAN')],
              'CL': tsv[tsv.rfind('/') + 1:tsv.rfind('.tsv')].strip('_')}[mode]

        print(tf, end=' ')

        # if not tf.startswith('ANDR'): break###

        with open(tsv) as file:

            for line in file:

                if not count:

                    title = line[:-1].split('\t')
                    count += 1

                else:

                    count += 1

                    a = line[:-1].split('\t')

                    snp = {title[i]: a[i] for i in range(len(a))}

                    if nagg == 'n' and float(snp['n_aggregated']) <= 1.0:
                        continue

                    try:
                        ref, alt = float(snp['%sref' % pv]), float(snp['%salt' % pv])
                    except:
                        ref, alt = float('NaN'), float('NaN')

                    try:
                        floatcallers = float(snp['n_peak_callers'])
                    except:
                        floatcallers = float(snp['m_callers'])

                    peakscond = {'>': floatcallers > 0,
                                 '==': floatcallers == 0,
                                 '>=': floatcallers >= 0}[peaks]

                    try:
                        s = int(snp['ID'][2:])

                        # if s==7873784: print(snp)

                        if min(ref, alt) < fdr and peakscond:

                            snps[s] = {'ref': set(), 'alt': set()}

                            if ref < fdr:
                                snps[s]['ref'].add(tf)
                            elif alt < fdr:
                                snps[s]['alt'].add(tf)
                    except:
                        continue
    return snps


def Make_snpphtf(mode, release_path, pval, pvallim, snpph_path, snpphtf_path):
    snps = Get_ASB(mode.upper(), release_path, pval, pvallim)

    print()
    print('Total %s ASB RS ids:' % mode.upper(), len(snps))

    outp2 = open(snpphtf_path, 'w')

    bnames = ['grasp', 'ebi', 'clinvar', 'phewas', 'finemapping']
    title = ['RSID', '#ref%ss' % mode, '#alt%ss' % mode,
             '#all', '#allbutgrasp', '#allsum', '#allsumbutgrasp'] + \
            ['#' + x for x in bnames] + ['ref%ss' % mode, 'alt%ss' % mode] + bnames

    outp2.write('\t'.join(title) + '\n')

    with open(snpph_path) as ffile:
        for line in ffile:
            if line.startswith('RSID'):
                tit = line.strip('\n').split('\t')
                titlen = len(tit)
                continue

            a = line.strip('\n').split('\t')
            a = {tit[x]: a[x] for x in range(titlen)}
            s = int(a['RSID'][2:])

            if s in snps:
                refalt = snps.pop(s)

                a['#ref%ss' % mode] = len(refalt['ref'])
                a['#alt%ss' % mode] = len(refalt['alt'])

                a['ref%ss' % mode] = ';'.join(sorted([y for y in refalt['ref']]))
                a['alt%ss' % mode] = ';'.join(sorted([y for y in refalt['alt']]))

                outp2.write('\t'.join([str(a[x]) for x in title]) + '\n')

    for s in snps:

        res = {}
        res['RSID'] = 'rs' + str(s)

        res['#all'] = 0
        res['#allbutgrasp'] = 0
        res['#allsum'] = 0
        res['#allsumbutgrasp'] = 0

        for bn in bnames:
            res['#' + bn] = 0
            res[bn] = ''

        res['#ref%ss' % mode] = len(snps[s]['ref'])
        res['#alt%ss' % mode] = len(snps[s]['alt'])

        res['ref%ss' % mode] = ';'.join(sorted([y for y in snps[s]['ref']]))
        res['alt%ss' % mode] = ';'.join(sorted([y for y in snps[s]['alt']]))

        outp2.write('\t'.join([str(res[x]) for x in title]) + '\n')

    outp2.close()

    print('{} is successfully created'.format(snpphtf_path))


def Concat_TF_CL(tffile, clfile, snpphtfcl_path):
    bnames = ['grasp', 'ebi', 'clinvar', 'phewas', 'finemapping']
    titlecl = ['RSID', '#refcls', '#altcls', '#all', '#allbutgrasp', '#allsum', '#allsumbutgrasp'] + \
              ['#' + x for x in bnames] + ['refcls', 'altcls'] + bnames
    titletf = ['RSID', '#reftfs', '#alttfs', '#all', '#allbutgrasp', '#allsum', '#allsumbutgrasp'] + \
              ['#' + x for x in bnames] + ['reftfs', 'alttfs'] + bnames

    title = ['RSID', '#reftfs', '#alttfs', '#refcls', '#altcls', '#all', '#allbutgrasp', '#allsum', '#allsumbutgrasp'] + \
            ['#' + x for x in bnames] + ['reftfs', 'alttfs'] + ['refcls', 'altcls'] + bnames

    snps = {}

    with open(tffile) as ffile:
        for line in ffile:
            if line.startswith('RSID'):
                tit = line.strip('\n').split('\t')
                titlen = len(tit)
                continue

            a = line.strip('\n').split('\t')
            a = {tit[x]: a[x] for x in range(titlen)}
            s = int(a['RSID'][2:])

            snps[s] = a

    outp = open(snpphtfcl_path, 'w')
    outp.write('\t'.join(title) + '\n')

    with open(clfile) as ffile:
        for line in ffile:
            if line.startswith('RSID'):
                tit = line.strip('\n').split('\t')
                titlen = len(tit)
                continue

            a = line.strip('\n').split('\t')
            a = {tit[x]: a[x] for x in range(titlen)}
            s = int(a['RSID'][2:])

            if s in snps:

                res = snps.pop(s)

                assert res['#allsum'] == a['#allsum']

                res['#refcls'] = a['#refcls']
                res['#altcls'] = a['#altcls']

                res['refcls'] = a['refcls']
                res['altcls'] = a['altcls']

                outp.write('\t'.join([str(res[x]) for x in title]) + '\n')

            else:

                res = a

                res['#reftfs'] = 0
                res['#alttfs'] = 0

                res['reftfs'] = ''
                res['alttfs'] = ''

                outp.write('\t'.join([str(res[x]) for x in title]) + '\n')

    for s in snps:
        res = snps[s]

        res['#refcls'] = 0
        res['#altcls'] = 0

        res['refcls'] = ''
        res['altcls'] = ''

        outp.write('\t'.join([str(res[x]) for x in title]) + '\n')

    outp.close()
    print('{} is successfully created'.format(snpphtfcl_path))


def Get_ASB_chr_pos(mode, release_path, pv='fdrp_bh_', fdr=0.05, peaks='>=', nagg=''):
    snps = {}

    files = glob.glob('{}/*.tsv'.format(release_path))

    for tsv in sorted(files):

        count = 0

        tf = {'TF': tsv[tsv.rfind('/') + 1:tsv.rfind('_HUMAN')],
              'CL': tsv[tsv.rfind('/') + 1:tsv.rfind('.tsv')].strip('_')}[mode]

        print(tf, end=' ')

        # if tf.startswith('A375__'): break###

        with open(tsv) as file:

            for line in file:

                if not count:

                    title = line[:-1].split('\t')
                    count += 1

                else:

                    count += 1

                    a = line[:-1].split('\t')

                    snp = {title[i]: a[i] for i in range(len(a))}

                    if nagg == 'n' and float(snp['n_aggregated']) <= 1.0:
                        continue

                    try:
                        ref, alt = float(snp['%sref' % pv]), float(snp['%salt' % pv])
                    except:
                        ref, alt = float('NaN'), float('NaN')

                    try:
                        floatcallers = float(snp['n_peak_callers'])
                    except:
                        floatcallers = float(snp['m_callers'])

                    peakscond = {'>': floatcallers > 0,
                                 '==': floatcallers == 0,
                                 '>=': floatcallers >= 0}[peaks]

                    try:
                        s = int(snp['ID'][2:])

                        # if s==7873784: print(snp)

                        if min(ref, alt) < fdr and peakscond:

                            if s not in snps:
                                snps[s] = (snp['#chr'], snp['pos'])

                    except:
                        continue
    return snps


def Add_chr_pos(release_TF_path, release_CL_path, pval, pvallim, tfclfile, output_path):
    snps = Get_ASB_chr_pos('TF', release_TF_path, pval, pvallim)
    snps.update(Get_ASB_chr_pos('CL', release_CL_path, pval, pvallim))

    bnames = ['grasp', 'ebi', 'clinvar', 'phewas', 'finemapping']
    title = ['RSID', 'chr', 'pos', '#reftfs', '#alttfs', '#refcls', '#altcls', '#all', '#allbutgrasp', '#allsum',
             '#allsumbutgrasp'] + \
            ['#' + x for x in bnames] + ['reftfs', 'alttfs'] + ['refcls', 'altcls'] + bnames

    outp = open(output_path, 'w')
    outp.write('\t'.join(title) + '\n')

    with open(tfclfile) as ffile:
        for line in ffile:
            if line.startswith('RSID'):
                tit = line.strip('\n').split('\t')
                titlen = len(tit)
                continue

            a = line.strip('\n').split('\t')
            a = {tit[x]: a[x] for x in range(titlen)}
            s = int(a['RSID'][2:])

            a['chr'] = snps[s][0]
            a['pos'] = snps[s][1]

            outp.write('\t'.join([str(a[x]) for x in title]) + '\n')

    outp.close()

    print()
    print('{} is successfully created'.format(output_path))


def Add_eQTL(qtlfiles, transqtl, tfclfile, output_path):
    title = ['RSID', 'chr', 'pos', '#reftfs', '#alttfs', '#refcls', '#altcls', '#QTL', '#QTLg', '#all', '#allbutgrasp',
             '#allsum', '#allsumbutgrasp'] + \
            ['#' + x for x in phenotype_db_names] + ['reftfs', 'alttfs'] + ['refcls', 'altcls'] + phenotype_db_names + [
                'QTL', 'QTLg']

    outp = open(output_path, 'w')
    outp.write('\t'.join(title) + '\n')

    snps = {}

    with open(tfclfile) as ffile:
        for line in ffile:
            if line.startswith('RSID'):
                tit = line.strip('\n').split('\t')
                titlen = len(tit)
                continue

            a = line.strip('\n').split('\t')
            a = {tit[x]: a[x] for x in range(titlen)}

            snps[a['chr'] + '_' + a['pos']] = [set(), set()]

    print('Number of cis-eQTL files:', len(qtlfiles))

    for qtlfile in qtlfiles:
        with open(qtlfile) as qfile:

            tis = qtlfile[qtlfile.rfind('/') + 1:qtlfile.find('.v8.')]
            print(tis, end=' ')

            for line in qfile:

                if line.startswith('variant_id'):
                    tit = line.strip('\n').split('\t')
                    titlen = len(tit)
                    continue

                a = line.strip('\n').split('\t')
                a = {tit[x]: a[x] for x in range(titlen)}

                chrpos = '_'.join(a['variant_id'].split('_')[:2])

                if chrpos in snps:
                    snps[chrpos][0].add(tis)
                    snps[chrpos][1].add(a['gene_id'])

    with open(transqtl) as trfile:
        for line in trfile:
            if line.startswith('tissue_id'):
                tit = line.strip('\n').split('\t')
                titlen = len(tit)
                continue

            a = line.strip('\n').split('\t')
            a = {tit[x]: a[x] for x in range(titlen)}

            chrpos = '_'.join(a['variant_id'].split('_')[:2])
            tis = a['tissue_id']
            gen = a['gene_id']

            if chrpos in snps:
                snps[chrpos][0].add(tis)
                snps[chrpos][1].add(gen)

    with open(tfclfile) as ffile:
        for line in ffile:
            if line.startswith('RSID'):
                tit = line.strip('\n').split('\t')
                titlen = len(tit)
                continue

            a = line.strip('\n').split('\t')
            a = {tit[x]: a[x] for x in range(titlen)}

            chrpos = a['chr'] + '_' + a['pos']

            a['#QTL'] = len(snps[chrpos][0])
            a['QTL'] = ';'.join(sorted(list(snps[chrpos][0])))

            a['#QTLg'] = len(snps[chrpos][1])
            a['QTLg'] = ';'.join(sorted(list(snps[chrpos][1])))

            outp.write('\t'.join([str(a[x]) for x in title]) + '\n')

    outp.close()
    print()
    print('{} is successfully created'.format(output_path))


def main(phenotypes_dir=None):
    if phenotypes_dir is None:
        phenotypes_dir = '/home/abramov/phenotypes'
    release_TF_path = get_result_dir_path('TF')
    release_CL_path = get_result_dir_path('CL')
    outp_path = os.path.join(get_release_stats_path(), 'phenotypes_stats.tsv')
    files = {'GRASP': os.path.join(phenotypes_dir, 'pheno', 'grasp_pheno.tsv'),
             'EBI': os.path.join(phenotypes_dir, 'pheno', 'gwas_catalog.tsv'),
             'ClinVar': os.path.join(phenotypes_dir, 'pheno', 'variant_summary.txt'),
             'FineMapping': os.path.join(phenotypes_dir, 'pheno', 'finemapping.csv'),
             'PheWas': os.path.join(phenotypes_dir, 'pheno', 'phewas-catalog.csv'),
             }
    pheno_tmp_dir = os.path.join(get_release_stats_path(), 'pheno_tmp')
    check_and_create_dir(pheno_tmp_dir)
    qtlfiles = glob.glob(os.path.join(phenotypes_dir, 'eqtl', 'signif', '*.txt'))
    transqtl = os.path.join(phenotypes_dir, 'eqtl', 'GTEx_Analysis_v8_trans_eGenes_fdr05.txt')

    pval = 'fdrp_bh_'
    pvallim = {'fdrp_bh_': 0.05, 'logitp_': 0.0005}[pval]

    snpph_path = os.path.join(pheno_tmp_dir, 'snpph.tsv')
    snpphtf_path = os.path.join(pheno_tmp_dir, 'snpphtfASB.tsv')
    snpphcl_path = os.path.join(pheno_tmp_dir, 'snpphclASB.tsv')
    snpphtfcl_path = os.path.join(pheno_tmp_dir, 'snpphtfclASB.tsv')
    snpphtfcl_chrpos_path = os.path.join(pheno_tmp_dir, 'snpphtfclASBchrpos.tsv')

    annotate_main(files, snpph_path)

    Make_snpphtf('tf', release_TF_path, pval, pvallim, snpph_path, snpphtf_path)
    Make_snpphtf('cl', release_CL_path, pval, pvallim, snpph_path, snpphcl_path)

    Concat_TF_CL(snpphtf_path, snpphcl_path, snpphtfcl_path)

    Add_chr_pos(release_TF_path, release_CL_path, pval, pvallim, snpphtfcl_path, snpphtfcl_chrpos_path)

    Add_eQTL(qtlfiles, transqtl, snpphtfcl_chrpos_path, outp_path)
    shutil.rmtree(pheno_tmp_dir)
