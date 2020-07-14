import sys
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt, ticker
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns

sns.set(font_scale=1.2, style="ticks", font="lato", palette=('#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2',
                                                             '#D55E00', '#CC79A7'))
# sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})
plt.rcParams['font.weight'] = "medium"
plt.rcParams['axes.labelweight'] = 'medium'
plt.rcParams['figure.titleweight'] = 'medium'
plt.rcParams['axes.titleweight'] = 'medium'
plt.rcParams['figure.figsize'] = 14, 6
plt.rcParams["legend.framealpha"] = 1
plt.rcParams['axes.xmargin'] = 0
plt.rcParams['axes.ymargin'] = 0
plt.rcParams["legend.framealpha"] = 1
plt.rcParams['figure.dpi'] = 432/6

get_color = {'A': '#0074FF', 'T': '#7900C8', 'G': '#FF4500', 'C': '#FFA500'}

for nuc in get_color:
    fig, ax = plt.subplots(figsize=(6, 0.6))
    fig.subplots_adjust(bottom=0.5)
    plt.tight_layout(w_pad=0, h_pad=0)
    top_pad = 0.02
    bottom_pad = 0.01
    text_h = 0.45
    left_pad = 0.023
    right_pad = 0.03
    ax.set_position([left_pad, text_h + bottom_pad, 1 - left_pad - right_pad, 1 - text_h - top_pad - bottom_pad])

    colors = ['white', get_color[nuc]]
    cmap = LinearSegmentedColormap.from_list("mycmap", colors)
    norm = mpl.colors.Normalize(vmin=0, vmax=20)
    cb = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                   norm=norm,
                                   orientation='horizontal')
    # ax.set_xlabel('-log10 FDR P-value of ASB')

    plt.savefig(os.path.expanduser('~/scale_{}.svg'.format(nuc)))
    plt.close(fig)
