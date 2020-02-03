import json
import os
import string


def remove_punctuation_bad(x):
    table = str.maketrans({key: "_" for key in string.punctuation if key not in '-'})
    return x.translate(table).replace(" ", "_")


def remove_punctuation(x):
    table = str.maketrans({key: "_" for key in string.punctuation if key not in '-+'})
    return x.translate(table).replace(" ", "_")


with open(os.path.expanduser('~/PARAMETERS/CELL_LINES.json')) as f1, open(
        os.path.expanduser('~/PARAMETERS/CELL_LINES_old.json')) as f2:
    d = json.loads(f1.readline())
    old_d = json.loads(f2.readline())

new_d = {}
new_converter = {}
convertPluses = {"CD8+ T-cells", "CD133+ umbilical cord blood cells",
                 "CD34+CD133+ hematopoietic progenitors", "CD4+ CD25- T-cells",
                 "CD4+ CD25+ CD45RA+ T-cells", "CD4+ T-cells", "CD34+ bone marrow cell-derived erythroid progenitors",
                 "CD34+ hematopoietic stem cells", "CD34+ hematopoietic stem cells-derived proerythroblasts",
                 "CD34+ hematopoietic stem progenitor cells", "CD34+ stem cells-derived erythroblasts",
                 "CD34+CD133+ hematopoietic progenitors", "fetal CD34+ hematopoietic stem progenitor cells",
                 "thymus-derived CD34+ cells", "CD14+ monocytes"}
d_items = d.items()
for element in convertPluses:
    new_converter[remove_punctuation_bad(element)] = remove_punctuation(element)
for key_old, value_old in old_d.items():
    for key, value in d_items:
        if value == value_old:
            for convCL in new_converter:
                if convCL != key.split("!")[0]:
                    continue
                print("Found ya {}".format(key))
                new_key = new_converter[key.split("!")[0]] + "!" + key.split("!")[1]
                new_d[key_old] = new_key
                #  change name in json too
                val = d[key]
                d[new_key] = val
                del d[key]
            else:
                new_d[key_old] = key
            break
    else:
        print('No match for {} ({})'.format(key_old, value_old))

# with open(os.path.expanduser('~/PARAMETERS/CELL_LINES.json'), "w") as f1:
#     json.dump(d, f1)
# path = os.path.expanduser('~/PloidyForRelease/CAIC/')
# for file in os.listdir(path):
#     print(path + file)
#     if file.replace('_ploidy.tsv', '') in new_d.values():
#         continue
#     os.rename(path + file, path + new_d[file.replace('_ploidy.tsv', '')] + '_ploidy.tsv')
#
# path = os.path.expanduser('~/PloidyForRelease/SQRT/')
# for file in os.listdir(path):
#     print(path + file)
#     if file.replace('_ploidy.tsv', '') in new_d.values():
#         continue
#     os.rename(path + file, path + new_d[file.replace('_ploidy.tsv', '')] + '_ploidy.tsv')
#
# path = os.path.expanduser('~/PloidyForRelease/merged_vcfs/')
# for file in os.listdir(path):
#     print(path + file)
#     if file.replace('.tsv', '') in new_d.values():
#         continue
#     os.rename(path + file, path + new_d[file.replace('.tsv', '')] + '.tsv')
