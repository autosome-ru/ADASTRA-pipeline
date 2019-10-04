import sys


table = open(sys.argv[1], 'r')
ape = open(sys.argv[2], 'r')
out = open(sys.argv[3], 'w')


#out.write('\t'.join(['#chr', 'pos', 'ID', 'ref', 'alt', 'm_callers', 'm_ploidy', 'm_q', 'm_dipq', 'm_segc', 'm_datasets', 'm_hpref', 'm_hpalt', 'm_fpref', 'm_fpalt', 'm_stpref', 'm_stpalt', 'prefectos_p1', 'perfectos_p2', 'perfectos_fc', 'motif_pos', 'orientation'])+'\n')

ape.readline()
ape_line = ape.readline()
table_line = table.readline()	
while ape_line and table_line:
	ape_line = ape_line.split()
	table_line = table_line.strip().split()
	if table_line[0][0] == "#":
		out.write("\t".join(table_line + ['perfectos_p1', 'perfectos_p2', 'perfectos_fc', 'motif_pos', 'orientation']) + "\n")
		table_line = table.readline().split()
	ape_name = ape_line[0].split('_')
	ape_chr = ape_name[0]
	ape_pos = ape_name[1]
	
	chr = table_line[0]
	pos = table_line[1]
	
	assert ape_chr == chr
	assert ape_pos == pos
	
	orientation = ape_line[3]

	if ape_line[2] == ape_line[5] and ape_line[3] == ape_line[6]:
		x = int(ape_line[2])
		Ln = len(ape_line[4])
		if orientation == 'direct':
			mpos = 1-x
		else:
			mpos = Ln+x-1
	else:
		mpos = '.'

	fc = float(ape_line[-1])
	
	p2 = round(float(ape_line[-2]), 50)
	p1 = round(float(ape_line[-3]), 50)
	out.write('\t'.join(['\t'.join(table_line), '\t'.join(map(str, [p1, p2, fc])), str(mpos), orientation])+'\n')
	
	ape_line = ape.readline()
	table_line = table.readline()

