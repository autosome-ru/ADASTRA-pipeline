from matplotlib import pyplot as plt
import matplotlib
from scipy import stats
import math
import sys
import os
import numpy as np
from PIL import Image
from pylab import *

def get_f(tup):
	strg = tup[0]
	f = tup[1] 
	if strg != '.':
		return f(strg)
	return f('0')

def signum(x):
	if x < 0:
		return -1
	elif x == 0:
		return 0
	else:
		return 1

def get_params(line):
	line = line.split()
	chr = line[0]
	pos = int(line[1])
	ID = line[2]
	ref = line[3]
	alt = line[4]
	(Q_exp, Q_ctrl) = map(float, line[5:7])
	(GQ_exp, GQ_ctrl, ref_c_exp, alt_c_exp, ref_c_ctrl, alt_c_ctrl, in_exp, in_ctrl, in_macs, in_sissrs, in_cpics, in_gem) = map(get_f, zip(line[7:19], [int]*12))
	(p_bin, p_fisher, pv_1, pv_2, fc) = map(get_f, zip(line[19:24], [float]*5))
	if line[24] != '.':
		(x, y) = map(int, line[24].split(':'))
		p_exp = x/(x+y)
	else: p_exp = 0
	if line[25] != '.':
		(x, y) = map(int, line[25].split(':'))
		p_ctrl = x/(x+y)
	else: p_ctrl = 0
	callers = in_macs + in_sissrs + in_cpics + in_gem
	if line[26] != '.':
		mpos = int(line[26])
	else: mpos = 0
	alignment = line[27]
	return chr, pos, ID, ref, alt, Q_exp, Q_ctrl, GQ_exp, GQ_ctrl, ref_c_exp, alt_c_exp, ref_c_ctrl, alt_c_ctrl, in_exp, in_ctrl, in_macs, in_sissrs, in_cpics, in_gem, callers, p_bin, p_fisher, pv_1, pv_2, fc, p_exp, p_ctrl, mpos, alignment

def comp(x):
	if x == 'A':
		return 'T'
	if x == 'T':
		return 'A'
	if x == 'G':
		return 'C'
	if x == 'C':
		return 'G'

def get_mutation_number(ref, alt):
	if ref == 'A' and alt == 'T':
		return (0,1)
	elif ref == 'T' and alt == 'A':
		return (0,-1)
	elif ref == 'A' and alt == 'G':
		return (1,1)
	elif ref == 'G' and alt == 'A':
		return (1,-1)
	elif ref == 'A' and alt == 'C':
		return (2,1)
	elif ref == 'C' and alt == 'A':
		return (2,-1)
	elif ref == 'T' and alt == 'G':
		return (3,1)
	elif ref == 'G' and alt == 'T':
		return (3,-1)
	elif ref == 'T' and alt == 'C':
		return (4,1)
	elif ref == 'C' and alt == 'T':
		return (4,-1)
	elif ref == 'G' and alt == 'C':
		return (5,1)
	elif ref == 'C' and alt == 'G':
		return (5,-1)

def get_position(mut_n, mpos, interval=10):
	return mpos*interval + mut_n + 1

def get_color(n):
	return 'C'+str(n%10)
def get_minor_pos(minor, mpos, interval=10):
	if minor == 'A':
		mn = 0
	elif minor == 'T':
		mn = 1
	elif minor == 'G':
		mn = 2
	elif minor == 'C':
		mn = 3
	return mpos*interval + mn*2

def get_text(mut_n):
	if mut_n == 0:
		return 'T/A'
	if mut_n == 1:
		return 'G/A'
	if mut_n == 2:
		return 'C/A'
	if mut_n == 3:
		return 'G/T'
	if mut_n == 4:
		return 'C/T'
	if mut_n == 5:
		return 'C/G'

def get_letter_color(L):
	if L == 'A': return 'DodgerBlue'
	if L == 'T': return 'Blueviolet'
	if L == 'G': return 'OrangeRed'
	if L == 'C': return 'Goldenrod'

def get_text_minor(mn):
	if mn == 0:
		return 'A'
	elif mn == 1:
		return 'T'
	elif mn == 2:
		return 'G'
	elif mn == 3:
		return 'C'

def get_text_split(mut_n, up):
	if up:
		return get_text(mut_n)[0]
	else:
		return get_text(mut_n)[-1]

def move_figure(f, x, y):
    """Move figure's upper left corner to pixel (x, y)"""
    backend = matplotlib.get_backend()
    if backend == 'TkAgg':
        f.canvas.manager.window.wm_geometry("+%d+%d" % (x, y))
    elif backend == 'WXAgg':
        f.canvas.manager.window.SetPosition((x, y))
    else:
        # This works for QT and GTK
        # You can also use window.setGeometry
        f.canvas.manager.window.move(x, y)

if __name__ == "__main__":
	pv_tr = 0.0005
	Q_tr = 0
	t_tr = 0
	count_tr = 0
	Q_max = 8
	trashy = 0
	pshy = 0
	Common = dict()
	for file in os.listdir("/home/sashok/Documents/ASB/CTCF_new_datasets/ALIGNS"):
		if file.endswith(".txt"):
			table = os.path.join("/home/sashok/Documents/ASB/CTCF_new_datasets/ALIGNS", file)
			print(table)
			table = open(table, 'r')
			for line in table:
				if line[0] == '#':
					continue
				(chr, pos, ID, ref, alt, Q_exp, Q_ctrl, GQ_exp, GQ_ctrl, ref_c_exp, alt_c_exp, ref_c_ctrl, alt_c_ctrl, in_exp, in_ctrl, in_macs, in_sissrs, in_cpics, in_gem, callers, p_bin, p_fisher, pv_1, pv_2, fc, p_exp, p_ctrl, mpos, alignment) = get_params(line)
				if callers < 1 or max(pv_1, pv_2) > pv_tr: continue
				if not in_exp: continue
				if Q_exp < Q_tr: continue
				if mpos==0: continue
				if alignment == 'revcomp':
					ref = comp(ref)
					alt = comp(alt)
					mpos += 1
				#if get_mutation_number(ref, alt)[1] == -1:
				#	(ref, alt) = (alt, ref)
				key = (chr, pos, ID, ref, alt, mpos, alignment)
				if p_exp == 0: continue
				if alt_c_exp < ref_c_exp:
					if ref_c_exp/(ref_c_exp+alt_c_exp) < 0.5 or ref_c_exp/(ref_c_exp+alt_c_exp) > 1: trashy += 1
					if 1 - p_exp < 0.5 or 1 - p_exp > 1: pshy += 1
					padd = math.log(ref_c_exp/(ref_c_exp+alt_c_exp)/(1 - p_exp), 2)
					if padd < 0: continue
					padd *= -1
					#padd = ref_c_exp/(ref_c_exp+alt_c_exp)/(1 - p_exp)
					#padd = 1/padd
				else:
					if alt_c_exp/(ref_c_exp+alt_c_exp) < 0.5 or alt_c_exp/(ref_c_exp+alt_c_exp) > 1: trashy += 1
					if 1 - p_exp < 0.5 or 1 - p_exp > 1: pshy += 1
					padd = math.log(alt_c_exp/(ref_c_exp+alt_c_exp)/(1 - p_exp), 2)
					if padd < 0: continue
					#padd = alt_c_exp/(ref_c_exp+alt_c_exp)/(1 - p_exp)
				#assert(padd<=2 and padd >=0.5)
				assert(padd>=-1 and padd <= 1)
				if Common.get(key, None):
					Common[key][p_bin] = padd#alt_c_exp/(ref_c_exp+alt_c_exp)
				else:
					Common[key] = dict()
					Common[key][p_bin] = padd#alt_c_exp/(ref_c_exp+alt_c_exp)
	sorted_Common = sorted(Common.items(), key=lambda x: len(x[1].keys()), reverse=True)
	out = open('/home/sashok/Documents/ASB/CTCF_stats/CTCF_common.txt', 'w')
	lengths = []
	colors = []
	markers = []
	fcs = []
	x = []
	y = []
	z=[]
	reported = 0
	cout = 0
	filt = 0
	interesting = []
	fisher = False
	log_pv_tr = 0
	for (key, pvs) in sorted_Common:
		list_pvs = []
		for pv in pvs.keys():
			list_pvs.append(pv)
			#colors.append((max(1-pvs[pv]/Q_max, 0), 0.3+0.4*(max(0, min(1, pvs[pv]/Q_max-1))), 0.8*min(pvs[pv]/Q_max, 1)))
		(chr, pos, ID, ref, alt, mpos, alignment) = key
		lengths.append(len(list_pvs))
		mean = np.mean(list_pvs)
		if len(list_pvs) >= count_tr:
			
			positive = 0
			negative = 0

			plus = []
			minus = []
			cout+=1
			AFs = []
			for pv in list_pvs:
				#x.append(pv)
				if pv > 0:
					plus.append(min(10**(-pv),1))
					positive += 1
				elif pv < 0:
					minus.append(min(10**(pv),1))
					negative += 1
				af = Common[key][pv]
				AFs.append(af)
			
			if plus != []:
				pluscomb = min(stats.combine_pvalues(plus)[1]*200000,1)
			else:
				pluscomb = 1
			if minus != []:
				minuscomb = min(stats.combine_pvalues(minus)[1]*200000,1)
			else:
				minuscomb = 1
			#print(pluscomb,minuscomb)
			#sign = signum(-pluscomb+minuscomb)
			#print(pluscomb,minuscomb,sign)
			if pluscomb <= 0 or minuscomb <= 0:
				reported += 1
				continue
			m = -1*math.log10(pluscomb/minuscomb)
			#if np.mean(AFs) <= 0:
			#	reported += 1				
			#	continue
			#AFm = math.log(np.mean(AFs), 2)
			AFm = np.mean(AFs)
			if abs(m) < log_pv_tr: continue
			if abs(m)>0:
				if abs(m) >100:
					interesting.append((mpos, chr, pos, ID, ref, alt, alignment))
				filt+=1
				#*get_mutation_number(ref, alt)[1])
				#y.append(max(AFm, 1-AFm))
				#y.append(fc)
				#z.append(key)
				if fisher:
					if abs(m) > 100: m *= 100/abs(m)
					y.append(abs(m))
					if m < 0:
						let = ref
						minor = alt
					else:
						let = alt
						minor = ref
				else:
					y.append(abs(AFm))		
					if AFm < 0:
						let = ref
						minor = alt
					else:
						let = alt
						minor = ref
				x.append(get_minor_pos(minor, mpos))
				colors.append(get_letter_color(let))
				#markers.append('$'+let+'$')
				

			fcs.append(fc)
	out=open("CTCF-fisher.txt","w")
	for i in range(len(z)):
		out.write("\t".join(map(str,z[i]))+ "\t" + str(x[i]) +"\n")
	print(cout,filt)
	#plt.hist(lengths)
	#plt.show()
	

	fig, ax = plt.subplots(figsize=(18.75,10.875))
	

	l, b, w, h = ax.get_position().bounds
	print(l,b,w,h)
	ax.set_position([l, b+h/9, w, h])
	

	inv = ax.transData.inverted()
	fwd = ax.transData

	print(fwd.transform([l, b])-fwd.transform([0, 0]))
	
	size = fwd.transform([l+w, b+h/9])-fwd.transform([l, b])
	xs = int(round(17.85/20*size[0]))
	ys = int(round(18/10*size[1]))

	print(size)
	xoff = fwd.transform([l, b])[0]-fwd.transform([0, 0])[0]
	yoff = fwd.transform([l, b])[1]-fwd.transform([0, 0])[1]

	im = Image.open("/home/sashok/Pictures/TVz02GxAPOc.jpg")
	im = im.resize((xs,ys), Image.ANTIALIAS)
	IM = fig.figimage(im, 216, 10)

	ax.set_xticklabels([])
	
	#for i in range(6):
	#	plt.text(10*(i+1)+3, 100, get_text(i), color='grey', size='large', weight='bold', ha='center')
	#	plt.text(10*(i+1)+3, 0.15, get_text(i), color=get_color(i), size='large', weight='bold', ha='center')
	#	plt.text(10*(i+1)+3, 100, get_text_split(i, True), color=get_color(i), size='large', weight='bold', ha='center')
	#	plt.text(10*(i+1)+3, -100, get_text_split(i, False), color=get_color(i), size='large', weight='bold', ha='center')
	#for i in range(19):
		#plt.text(10*(i+1)+3, 150, str(i+1), color='grey', size='large', ha='center')
		#plt.text(10*(i+1)+3, 0.25, str(i+1), color='grey', size='large', ha='center')
	if fisher:
		plt.scatter(x, y, c=colors, s=15)
		plt.ylim(-20, 101)
	else:
		plt.scatter(x, y, c=colors, s=15)
		plt.ylim(-1.1, 1.1)
	#plt.text(mean+0.01, fc+0.01, ID, va="bottom", ha="left")
	plt.title( "perfectos_pv <= {}, log10(p_combined) >= {}".format(pv_tr, log_pv_tr))
	#plt.xlabel("position-mutation")
	if fisher:
		plt.ylabel("-log10(p)")
	else:
		plt.ylabel("log2(MajorAlleleFrequency/(1-p_estimated))")
	if fisher:
		h = -10
		hn = 75
	else:
		h = -0.25
		hn = -0.5
		plt.text(10*(1)+3, h*0.75, 'minors', color='grey', ha='center')
	for i in range(20):
		plt.axvline(x=10*i+8, linewidth=1, linestyle='--', color='grey')
		if i==19: continue
		plt.text(10*(i+1)+3, hn, str(i+1), color='grey', size='large', ha='center')
		for j in range(4):
			plt.axvline(x=10*i+10+2*j, linewidth=0.5, linestyle='--', color='grey', alpha=0.5)
			plt.text(10*i+10+2*j, h, get_text_minor(j), color=get_letter_color(get_text_minor(j)), ha='center')
	plt.axhline(y=0, linewidth=1, linestyle='--', color='grey')
	plt.axhline(y=-1, linewidth=1, linestyle='--', color='grey')
	if not fisher:
		#plt.axhline(y=0.5, linewidth=1, linestyle='--', color='grey')
		plt.axhline(y=1, linewidth=1, linestyle='--', color='grey')
	#plt.scatter(means, fcs, s=100, c=mean_colors)
	for i in interesting:
		print(i)
	print('reported {}, trashy {}, fishy {}'.format(reported, trashy, pshy))
	plt.show()
	
