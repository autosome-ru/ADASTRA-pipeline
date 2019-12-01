import os
import sys
import gzip

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths import make_black_list, parameters_path

with open(parameters_path + "Master-lines.tsv", 'r') as file:
    blacklist = make_black_list()
    controls = set()
    for line in file:
        line = line.strip().split('\t')
        if line[0] not in blacklist:
            al_path = "/home/abramov/Alignments/" + "EXP/" + line[1] + "/" + line[0] + "/" + line[6] + '.vcf.gz'
            chr11_path = '/home/abramov/Alignments_chr11/' + "EXP/" + line[1] + "/" + line[0] + "/" + line[
                6] + '.vcf.gz'

            a = os.path.isfile(al_path)
            c = os.path.isfile(chr11_path)
            if a and not c:
                print('No CHR11: {}'.format(chr11_path))
            elif c and not a:
                print('No ALIGN: {}'.format(al_path))
            elif not c and not a:
                print('No: {}'.format(al_path))
            else:
                try:
                    with gzip.open(chr11_path, 'rt') as chr11, gzip.open(al_path, 'at') as al:
                        for line in chr11:
                            if line[0] == '#':
                                continue
                            al.write(line)
                except:
                    print("EXCEPTION {}, {}".format(e.args[0], al_path))

        if len(line) >= 10 and line[10] not in blacklist:
            al_path = "/home/abramov/Alignments/" + "CTRL/" + line[10] + "/" + line[14] + '.vcf.gz'
            chr11_path = '/home/abramov/Alignments_chr11/' + "CTRL/" + line[10] + "/" + line[14] + '.vcf.gz'

            if al_path in controls:
                continue

            a = os.path.isfile(al_path)
            c = os.path.isfile(chr11_path)
            if a and not c:
                print('No CHR11: {}'.format(chr11_path))
            elif c and not a:
                print('No ALIGN: {}'.format(al_path))
            elif not c and not a:
                print('No: {}'.format(al_path))
            else:
                try:
                    with gzip.open(chr11_path, 'rt') as chr11, gzip.open(al_path, 'at') as al:
                        for line in chr11:
                            if line[0] == '#':
                                continue
                            al.write(line)
                except Exception as e:
                    print("EXCEPTION {}, {}".format(e.args[0], al_path))
            controls.add(al_path)

