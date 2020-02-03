import json
import os

with open(os.path.expanduser('~/PARAMETERS/CELL_LINES.json')) as f1, open(os.path.expanduser('~/PARAMETERS/CELL_LINES_old.json')) as f2:
    d = json.loads(f1.readline())
    old_d = json.loads(f2.readline())

new_d = {}

for key_old, value_old in old_d.items():
    if "+" in key_old:
        print(key_old)
    for key, value in d.items():

        if value == value_old:
            new_d[key_old] = key
            break
    else:
        print('No match for {} ({})'.format(key_old, value_old))

path = os.path.expanduser('~/PloidyForRelease/CAIC/')
for file in os.listdir(path):
    print(path + file)
    if file.replace('_ploidy.tsv', '') in new_d.values():
        continue
    os.rename(path + file, path + new_d[file.replace('_ploidy.tsv', '')] + '_ploidy.tsv')

path = os.path.expanduser('~/PloidyForRelease/SQRT/')
for file in os.listdir(path):
    print(path + file)
    if file.replace('_ploidy.tsv', '') in new_d.values():
        continue
    os.rename(path + file, path + new_d[file.replace('_ploidy.tsv', '')] + '_ploidy.tsv')

path = os.path.expanduser('~/PloidyForRelease/merged_vcfs/')
for file in os.listdir(path):
    print(path + file)
    if file.replace('.tsv', '') in new_d.values():
        continue
    os.rename(path + file, path + new_d[file.replace('.tsv', '')] + '.tsv')
