import os
from scripts.HELPERS.helpers import pack
from scripts.HELPERS.paths import get_release_stats_path

phenotype_db_names = ['grasp', 'ebi', 'clinvar', 'phewas', 'finemapping']

def parse_grasp(filepath):
    phenotypes = {}
    with open(filepath) as file:
        for line in file:
            a = line.strip('\n').split('\t')
            if a[1] not in phenotypes:
                phenotypes[a[1]] = set()
            phenotypes[a[1]].add(int(a[0]))

    return phenotypes


def parse_finemapping(filepath):
    phenotypes = {}
    with open(filepath) as file:
        for k, line in enumerate(file):
            if k == 0:
                continue
            a = line.strip('\n').split(',')
            if a[0] not in phenotypes:
                phenotypes[a[0]] = set()
            phenotypes[a[0]].add(int(a[2][2:]))

    return phenotypes


def parse_ebi(filepath):
    phenotypes = {}
    with open(filepath, 'r') as file:
        for k, line in enumerate(file):
            if k == 0:
                continue
            try:
                a = line.strip('\n').split('\t')
                if a[7] not in phenotypes:
                    phenotypes[a[7]] = set()
                phenotypes[a[7]].add(int(a[23]))
            except ValueError:
                continue
    return phenotypes


def parse_phewas(filepath):
    phenotypes = {}
    with open(filepath, 'r') as file:
        for k, line in enumerate(file):
            if k == 0:
                continue
            a = line[line.find('"rs'):].split('",')
            ph = a[1][1:]
            if ph not in phenotypes:
                phenotypes[ph] = set()
            phenotypes[ph].add(int(a[0][3:]))

    return phenotypes


def parse_clinvar(filepath):
    phenotypes = {}
    with open(filepath, 'r') as file:
        for k, line in enumerate(file):
            if k == 0:
                continue
            a = line.strip('\n').split('\t')
            if 'pathogenic' not in a[6].lower() and 'risk factor' not in a[6].lower():
                continue

            for ph in a[13].split(';'):
                if ph in ('not provided', 'not specified'):
                    continue
                if ph not in phenotypes:
                    phenotypes[ph] = set()
                phenotypes[ph].add(int(a[9]))
    return phenotypes


def main(files, path_to_output):
    phenotypes_for_db_list = [parse_grasp(files['GRASP']),
                              parse_ebi(files['EBI']),
                              parse_clinvar(files['ClinVar']),
                              parse_phewas(files['PheWas']),
                              parse_finemapping(files['FineMapping'])]

    phenotypes_ids_dict = {}
    ids_phenotypes_dict = {}
    phenotype_id = 1

    def remove_phen_name_punctuation(phenotype_name):
        return phenotype_name.lower().replace("'", '').replace('_', ' ')

    for db in phenotypes_for_db_list:
        for phenotype in db:
            if remove_phen_name_punctuation(phenotype) not in phenotypes_ids_dict:
                phenotypes_ids_dict[remove_phen_name_punctuation(phenotype)] = phenotype_id
                ids_phenotypes_dict[phenotype_id] = remove_phen_name_punctuation(phenotype)
                phenotype_id += 1

    all_phenotypes = {}

    for i in range(len(phenotypes_for_db_list)):
        for phenotype in phenotypes_for_db_list[i]:
            for s in phenotypes_for_db_list[i][phenotype]:
                if s not in all_phenotypes:
                    all_phenotypes[s] = {x: set() for x in phenotype_db_names}
                all_phenotypes[s][phenotype_db_names[i]].add(
                    phenotypes_ids_dict[remove_phen_name_punctuation(phenotype)])

    print('pheno sizes:', len(phenotypes_ids_dict), len(all_phenotypes))
    with open(path_to_output, 'w') as out:
        header = ['RSID', '#all', '#allbutgrasp', '#allsum', '#allsumbutgrasp'] + \
                ['#' + x for x in phenotype_db_names] + \
                phenotype_db_names

        out.write('\t'.join(header) + '\n')
        for s in all_phenotypes:
            abn, phenotypes_without_grasp = set(), set()
            for db in phenotype_db_names:
                if db != 'grasp':
                    phenotypes_without_grasp.update(all_phenotypes[s][db])
                abn.update(all_phenotypes[s][db])
            bb = [len(all_phenotypes[s][x]) for x in phenotype_db_names]
            bb = [sum(len(all_phenotypes[s][x]) for x in phenotype_db_names),
                  sum(len(all_phenotypes[s][x]) for x in phenotype_db_names[1:])] + bb
            bb = [len(abn), len(phenotypes_without_grasp)] + bb
            cc = [';'.join(sorted([ids_phenotypes_dict[y]
                                   for y in all_phenotypes[s][x]])) for x in phenotype_db_names]
            out.write(pack(['rs{}'.format(s), *bb, *cc]))
            for x in phenotype_db_names:
                all_phenotypes[s][x] = len(all_phenotypes[s][x])
    print('{} is successfully created'.format(path_to_output))
