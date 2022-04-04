import multiprocessing as mp
import os
import subprocess
import json

from scripts.HELPERS.helpers import get_merged_badmaps_dict_path, get_results_file
from scripts.HELPERS.paths import get_release_stats_path, get_ending, get_dir_by_stage
from scripts.HELPERS.paths_for_components import results_path


def read_badmaps(remade):
    with open(get_merged_badmaps_dict_path(remade=remade), 'r') as f:
        return json.load(f)


def create_filename_list_mixalime(badmap, exps, out):
    f_name = os.path.join(out, f'{badmap}.tables.txt')
    with open(f_name, 'w') as f:
        for exp in exps:
            path = get_results_file(exp, 'BAD')
            if os.path.exists(path):
                f.write(f'{path}\n')
            else:
                print(path)
    return f_name


def process_dataset(data):
    badmap_name, file_path = data
    print(f'Processing {badmap_name}')
    out_dir = get_dir_by_stage('p-value')
    pr = subprocess.run(['calc_pval', '-f', file_path, '-w', get_release_stats_path(),
                         '-O', out_dir, '-m', 'window'])


def main(remade=True, n_jobs=1):
    badmaps = read_badmaps(remade)
    out = os.path.join(results_path, 'pvalue_files')
    if not os.path.exists(out):
        os.mkdir(out)
    out_dir = get_dir_by_stage('p-value')
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    files_dict = {badmap: create_filename_list_mixalime(badmap, exps, out)
                  for badmap, exps in badmaps.items()}

    ctx = mp.get_context("forkserver")
    with ctx.Pool(n_jobs) as p:
        p.map(process_dataset, [(x, files_dict[x]) for x in badmaps])


if __name__ == '__main__':
    main()
