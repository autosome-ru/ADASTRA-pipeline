import sys

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.helpers import pack

motif_length = int(sys.argv[4])
#read sarus file and choose best
dict_of_snps = {}

with open(sys.argv[2], 'r') as sarus:
    allele = None
    current_snp_id = None
    for line in sarus:
        if line[0] == ">":
            # choose best
            allele = line[-4:-1]
            current_snp_id = line[1:-5]
            assert allele in ("ref", "alt")
            if allele == "ref":
                dict_of_snps[current_snp_id] = {"ref": [], "alt": []}
        else:
            assert allele in ("ref", "alt")
            line = line.split("\t")
            dict_of_snps[current_snp_id][allele].append({
                "p": float(line[0]),
                "orientation": line[2],
                "pos": int(line[1]) if line[2] == '-' else motif_length - 1 - int(line[1]),
            })

with open(sys.argv[1], 'r') as table, open(sys.argv[3], 'w') as out:
    for line in table:
        line = line.strip('\n').split('\t')
        if line[0][0] == '#':
            out.write(pack(line + ['motif_log_pref', 'motif_log_palt', 'fold_change', 'motif_pos', 'orientation']))
            continue
        ID = line[2]

        dict_of_snps[ID]['ref'] = sorted(dict_of_snps[ID]['ref'], key=lambda x: x['pos'])
        dict_of_snps[ID]['ref'] = sorted(dict_of_snps[ID]['ref'], key=lambda x: x['orientation'])

        dict_of_snps[ID]['alt'] = sorted(dict_of_snps[ID]['alt'], key=lambda x: x['pos'])
        dict_of_snps[ID]['alt'] = sorted(dict_of_snps[ID]['alt'], key=lambda x: x['orientation'])

        ref_best = max(enumerate(dict_of_snps[ID]['ref']), key=lambda x: x[1]['p'])
        alt_best = max(enumerate(dict_of_snps[ID]['alt']), key=lambda x: x[1]['p'])

        best_idx, _ = max((ref_best, alt_best), key=lambda x: x[1]['p'])

        assert len(dict_of_snps[ID]['ref']) == len(dict_of_snps[ID]['alt'])
        assert dict_of_snps[ID]['ref'][best_idx]['pos'] == dict_of_snps[ID]['alt'][best_idx]['pos']

        out.write(pack(line + [dict_of_snps[ID]['ref'][best_idx]['p'],
                               dict_of_snps[ID]['alt'][best_idx]['p'],
                               dict_of_snps[ID]['alt'][best_idx]['p'] - dict_of_snps[ID]['ref'][best_idx]['p'],
                               dict_of_snps[ID]['ref'][best_idx]['pos'],
                               dict_of_snps[ID]['ref'][best_idx]['orientation'],
                               ]))
