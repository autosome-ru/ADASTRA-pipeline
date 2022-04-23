import os
import subprocess
import json
from scripts.HELPERS.paths import get_release_stats_path


def main(dist):
    in_dir = get_release_stats_path()
    with open(os.path.join(in_dir, 'mixalime_params.json')) as f:
        mixalime_params = json.load(f)
    args = (item for pair in mixalime_params.items() for item in pair)
    pr = subprocess.run(['negbin_fit', '-O', in_dir, '-d', dist, *args],
                        check=True,
                        capture_output=True)
    err = pr.stderr
    if err:
        err = err.decode('utf-8')
        with open(os.path.join(in_dir, 'mixalime_log.txt'), 'w') as f:
            f.write(err)
