import sys
import os
import numpy as np
import pandas as pd
sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.FIGURES import style_config
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns

#for file_name in os.listdir('.'):
	#if not file_name.endswith('_union.tsv'):
	#	continue
	#cl = file_name[:file_name.find('_union.tsv')]
	## Import Data
	#df = pd.read_table(file_name, header=None)
	#df.columns = ['chr', 'pos', 'BAD', 'COSMIC']
	#df.loc[df.COSMIC >= 7, 'COSMIC'] = 7
	#df_counts = df.groupby(['BAD', 'COSMIC']).size().reset_index(name='counts')

if __name__ == '__main__':
	file_name = 'all_lines_union.tsv'
	cl = 'all_lines no filter'
	df_counts = pd.read_table(os.path.expanduser('~/Documents/ASB/Correlation/counts.tsv'))
	# Draw Stripplot
	fig, ax = plt.subplots(figsize=(16,10), dpi= 80)   
	ax.set_yticks(list(set(df_counts.BAD)), minor=False)
	ax.set_xticks(list(set(df_counts.COSMIC)), minor=False)
	ax.yaxis.grid(True, which='major')
	ax.xaxis.grid(True, which='major')
	ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

	sns.scatterplot(df_counts.COSMIC, df_counts.BAD, size=df_counts.counts**0.5, sizes=(10, 4000),
		hue=df_counts.COSMIC, palette=sns.color_palette('husl', len(set(df_counts.COSMIC))),
		legend=False)

	ax.set_axisbelow(True)

	# Decorations
	plt.title('Counts Plot - {}'.format(cl), fontsize=22)
	plt.savefig(os.path.expanduser('~/Documents/ASB/Correlation/Counts_plot_{}.png'.format(cl)))
	plt.show()
