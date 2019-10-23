import sys
import gzip

Nucleotides = {'A', 'T', 'G', 'C'}


def read_from_file(vcf, out):
    skipped = 0
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
                if GT != '0/1':
                    skipped += 1
                    continue
                out[(line[0], line[1])] = [R, A, NAME, REF, ALT, QUAL, GQ, GT]
    print('Skipped {} homozigous SNPs'.format(skipped))


def Write(in_key, num, ord_dictionary):
    ord_dictionary[in_key] = in_key + [num]


def less(A, B):
    if A[0] < B[0]:
        return True
    elif A[0] == B[0]:
        return A[1] < B[1]
    else:
        return False


def write_peak(in_key, peak_line, ord_dictionary, gem=False):
    chr = in_key[0]
    chr_p = "chr" + peak_line[0]
    pos = int(in_key[1])
    start = int(peak_line[1])
    end = int(peak_line[2])

    if gem:
        start = max(0, start - 200)
        end = end + 200

    if less((chr, pos), (chr_p, end)):
        if not less((chr, pos), (chr_p, start)):
            Write(in_key, '1', ord_dictionary)
        else:
            Write(in_key, '0', ord_dictionary)
        return True
    else:
        return False


def add_caller(caller, ord_dictionary, gem=False):
    peak_line = caller.readline()
    while peak_line and peak_line[0] == "#":
        peak_line = caller.readline()
    peak_line = peak_line.split()

    i = 0
    while i < len(ord_dictionary):
        while peak_line:
            if not write_peak(key, peak_line, ord_dictionary, gem):
                peak_line = caller.readline().split()
            else:
                in_line = input_file.readline().split()
        while in_line:
            Write(in_line, '0', output)
            in_line = input_file.readline().split()
    output.close()


def add_zeros(array):
    for key in array:
        array[key] = ord_dictionary[key] + ["0"]


# TODO FIX Annotate.py

if __name__ == "__main__":
    vcf = gzip.open(sys.argv[1], 'rt')

    exp = dict()
    read_from_file(vcf, exp)

    exp_keys = list(exp.keys())
    exp_keys = sorted(exp_keys, key=lambda x: x[1])
    exp_keys = sorted(exp_keys, key=lambda x: x[0])

    out = sys.argv[10]

    with_macs = sys.argv[6]
    with_sissrs = sys.argv[7]
    with_cpics = sys.argv[8]
    with_gem = sys.argv[9]

    if with_macs == "true":
        macs = open(sys.argv[2], "r")
        add_caller(macs, exp_keys)
    else:
        add_zeros(exp_keys)

    if with_sissrs == "true":
        sissrs = open(sys.argv[3], "r")
        add_caller(sissrs, exp_keys)
    else:
        add_zeros(exp_keys)

    if with_cpics == "true":
        cpics = open(sys.argv[4], "r")
        add_caller(cpics, exp_keys)
    else:
        add_zeros(exp_keys)

    if with_gem == "true":
        gem = open(sys.argv[5], "r")
        add_caller(gem, exp_keys, gem=True)
    else:
        add_zeros(exp_keys)
