import json
import os

with open(os.path.expanduser('~/PARAMETERS/CELL_LINES.json')) as f1, open(os.path.expanduser('~/PARAMETERS/CELL_LINES_old.json')) as f2:
    d = json.loads(f1.readline())
    old_d = json.loads(f2.readline())

new_d = {}
for key, value in d.items():
    for key_old, value_old in old_d.items():
        if value == value_old:
            new_d[key_old] = key
            break

path = '~/PloidyForRelease/CAIC/'
for file in os.listdir(os.path.expanduser(path)):
    os.rename(path + file, path + new_d[file.replace('_ploidy.tsv', '')] + '_ploidy.tsv')

path = '~/PloidyForRelease/SQRT/'
for file in os.listdir(os.path.expanduser(path)):
    os.rename(path + file, path + new_d[file.replace('_ploidy.tsv', '')] + '_ploidy.tsv')

path = '~/PloidyForRelease/merged_vcfs/'
for file in os.listdir(os.path.expanduser(path)):
    os.rename(path + file, path + new_d[file.replace('.tsv', '')] + '.tsv')
