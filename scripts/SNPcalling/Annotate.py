import sys
import gzip


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
                out[(line[0], line[1])] = (R, A, NAME, REF, ALT, QUAL, GQ, GT)
    print('Skipped {} homozigous SNPs'.format(skipped))


def Write(line, num, output):
    output.write('\t'.join(line) + '\t' + num + '\n')


def less(A, B):
    if A[0] < B[0]:
        return True
    elif A[0] == B[0]:
        return A[1] < B[1]
    else:
        return False


def write_peak(in_line, peak_line, output, gem=False):
    chr = in_line[0]
    chr_p = "chr" + peak_line[0]
    pos = int(in_line[1])
    start = int(peak_line[1])
    end = int(peak_line[2])

    if gem:
        start = max(0, start - 200)
        end = end + 200

    if less((chr, pos), (chr_p, end)):
        if not less((chr, pos), (chr_p, start)):
            Write(in_line, '1', output)
        else:
            Write(in_line, '0', output)
        return True
    else:
        return False


def add_caller(caller, infile, outfile, gem=False):
    input_file = open(infile, "r")
    output = open(outfile, "w")

    in_line = input_file.readline()

    while in_line and in_line[0] == '#':
        in_line = input_file.readline()
    in_line = in_line.split()

    peak_line = caller.readline()

    while peak_line and peak_line[0] == "#":
        peak_line = caller.readline()
    peak_line = peak_line.split()

    while in_line and peak_line:
        if not write_peak(in_line, peak_line, output, gem):
            peak_line = caller.readline().split()
        else:
            in_line = input_file.readline().split()
    while in_line:
        Write(in_line, '0', output)
        in_line = input_file.readline().split()
    input_file.close()
    output.close()


def add_zeros(name, infile, outfile):
    inp = open(infile, "r")
    output = open(outfile, "w")
    output.write("#No {} peaks".format(name))
    for line in inp:
        if line[0] == '#':
            output.write(line)
            continue
        output.write(line[:-1] + '\t' + '0' + '\n')
    inp.close()
    output.close()


Nucleotides = {'A', 'T', 'G', 'C'}
#TODO FIX Annotate.py
if __name__ == "__main__":
    vcf = gzip.open(sys.argv[1], 'rt')

    exp = dict()
    read_from_file(vcf, exp)

    exp_keys = list(exp.keys())
    exp_keys = sorted(exp_keys, key=lambda x: x[1])
    exp_keys = sorted(exp_keys, key=lambda x: x[0])

    out = sys.argv[10]
    
    withmacs = sys.argv[6]
    withsissrs = sys.argv[7]
    withcpics = sys.argv[8]
    withgem = sys.argv[9]
    
    if withmacs == "true":
        macs = open(sys.argv[2], "r")
        add_caller(macs, i, out + ".m.txt")
    else:
        add_zeros("macs", i, out + ".m.txt")
    
    if withsissrs == "true":
        sissrs = open(sys.argv[3], "r")
        add_caller(sissrs, out + ".m.txt", out + ".s.txt")
    else:
        add_zeros("sissrs", out + ".m.txt", out + ".s.txt")
    
    if withcpics == "true":
        cpics = open(sys.argv[4], "r")
        add_caller(cpics, out + ".s.txt", out + ".c.txt")
    else:
        add_zeros("cpics", out + ".s.txt", out + ".c.txt")
    
    if withgem == "true":
        gem = open(sys.argv[5], "r")
        add_caller(gem, out + ".c.txt", out, gem=True)
    else:
        add_zeros("gem", out + ".c.txt", out)
