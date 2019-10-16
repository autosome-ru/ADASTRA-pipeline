import sys
import gzip


def write(chr, pos, NAME, REF, ALT, R, A, Q, GQ, output):
    R = str(R)
    A = str(A)
    output.write(chr + '\t' + pos + '\t' + NAME + '\t' + REF + '\t' + ALT + '\t'
                 + R + '\t' + A + '\t' + Q + '\t' + GQ + '\n')


vcf = gzip.open(sys.argv[1], 'rt')
output = open(sys.argv[2], 'w')

exp = dict()

Nucleotides = {'A', 'T', 'G', 'C'}


def read_from_file(vcf, out):
    for line in vcf:
        if line[0] == '#':
            continue
        line = line.split()
        if len(line[3]) == 1 and len(line[4]) == 1:
            if line[3] in Nucleotides and line[4] in Nucleotides:
                Inf = line[-1].split(':')
                R = int(Inf[1].split(',')[0])
                A = int(Inf[1].split(',')[1])
                GQ = Inf[3]
                GT = Inf[0]
                NAME = line[2]
                REF = line[3]
                ALT = line[4]
                QUAL = line[5]
                out[(line[0], line[1])] = (R, A, NAME, REF, ALT, QUAL, GQ, GT)


read_from_file(vcf, exp)

skipped = 0

exp_keys = list(exp.keys())
exp_keys = sorted(exp_keys, key=lambda x: x[1])
exp_keys = sorted(exp_keys, key=lambda x: x[0])

for (chr, pos) in exp_keys:
    (R, A, NAME, REF, ALT, Q, GQ, GT) = exp[(chr, pos)]

    if GT != '0/1':
        skipped += 1
        continue
    else:
        write(chr, pos, NAME, REF, ALT, R, A, Q, GQ, output)

print('Skipped {} homozigous SNPs'.format(skipped))
