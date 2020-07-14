import sys
import numpy as np
from matplotlib import pyplot as plt

plt.tight_layout()
markers = ['o', (4, 0), (5, 0), (5, 1), (6, 1), 'D', 'x', '+', (5, 2), (6, 2)]

with open(sys.argv[1], 'r') as file:
	file_c = -1
	for path in file:
		file_c += 1
		path = path.strip()
		with open(path, 'r') as table:
			print(path,path.rfind('/'),path.find('_'))
			name = path[path.rfind('/'):path.find('_')]
			lines = dict()
			for line in table:
				lin = line.split()
				ref = int(lin[5])
				alt = int(lin[6])
				n = ref + alt
				rs = lin[2]
				cne = int(lin[7])
				if lin[9] == '0':
					cnv = 0
				else:
					cnv = int(lin[10])/int(lin[9])-1
				if cne != 1 or cnv != 1 or rs == '.':
					continue
				try:
					lines[n].append(ref)
				except KeyError:
					lines[n] = [ref]
			x = list(range(6, 51))
			y = []
			mean = 0
			for n in range(6, 51):
				f = 2*np.mean(lines[n])/n
				mean += f/45
				try:
					y.append(f)
				except KeyError:
					y.append(0)
			plt.scatter(x, y, label=name, s=40, marker=markers[file_c], color='C'+str(file_c))
			plt.hlines(mean, 1, 50, color='C'+str(file_c), linestyle='--', lw=1, label=name+' average')
			

x = list(range(6, 51))
y = [1 for i in x]
plt.plot(x, y, color='Black', label='Symmetrical')
plt.grid(True)
plt.title('Relative Ref bias in diploid regions')
plt.xlabel('Cover')
plt.ylabel('(Average Ref counts)/(cover/2)')
plt.legend()
plt.show()
