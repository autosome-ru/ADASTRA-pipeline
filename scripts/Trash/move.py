import shutil
import os
import sys

exp = os.path.expanduser('~/ALIGNMENTS/EXP')
ctrl = os.path.expanduser('~/ALIGNMENTS/CTRL')

out = os.path.expanduser(sys.argv[1])
if not os.path.isdir(out):
    print('Out is not a dir')
    exit(1)

for TF in os.listdir(exp):
    for exp in os.listdir(os.path.join(exp, TF)):
        if os.path.isdir(os.path.join(exp, TF, exp)):
            for file in os.path.join(exp, TF, exp):
                if file.endswith('.vcf.gz'):
                    break
            else:
                continue
            shutil.move(os.path.join(exp, TF, exp), os.path.join(out, exp))

for exp in os.listdir(ctrl):
    if os.path.isdir(os.path.join(ctrl, exp)):
        shutil.move(os.path.join(ctrl, exp), os.path.join(out, exp))
