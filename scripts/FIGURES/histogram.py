from scipy import special
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
from numpy.random import *
import sys

table = sys.argv[1]
name = table[table.rfind('/')+1:table.find('_')]
chrom = sys.argv[2]
N = int(sys.argv[3])

X = []
with open(table, 'r') as file:
	for line in file:
		if line[0] == '#':
			continue
		line = line.split()
		chr = line[0]
		rs = line[2]
		if rs == '.':
			continue
		cne = int(line[7])
		cnm = int(line[-2])
		cns = int(line[-1])
		if cnm == 0:
			cnv = 0
		else:
			cnv = cns/cnm - 1
		if chrom == 'all_dip_cnv':
			if cnv != 1:
				continue
			col = 'C0'
		elif chrom == 'all_dip_estimations':
			if cne != 1:
				continue
			col = 'C1'
		elif chrom == 'all_poly_cnv':
			if cnv <= 1:
				continue
			col = 'C0'
		elif chrom == 'all_poly_estimations':
			if cne <= 1:
				continue
			col = 'C1'
		ref = int(line[5])
		alt = int(line[6])
		x = ref
		n = ref+alt
		if n != N:
			continue
		X.append(x)

X = sorted(X)
print(X)
plt.hist(X,bins=N+1,range=(-0.5,N+0.5), color=col)
plt.grid()
plt.xlabel('ref counts')
plt.title('Hist, '+chrom+', n='+str(N)+', '+name)
plt.savefig('HIST_'+name+'_'+chrom+'_n='+str(N)+'.png', bbox_inches='tight')
plt.show()

