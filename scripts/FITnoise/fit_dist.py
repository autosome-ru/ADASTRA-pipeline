import os
import subprocess

from scripts.HELPERS.paths import get_release_stats_path


def main(model='NB_AS'):
    in_dir = get_release_stats_path()
    pr = subprocess.run(['negbin_fit', '-O', in_dir, '-m', model])
    err = pr.stderr.decode('utf-8')
    if err:
        print(err)
