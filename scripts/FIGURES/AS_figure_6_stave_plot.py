import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, ticker
import seaborn as sns
from PIL import Image

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.helpers import states

sns.set(font_scale=2.4, style="ticks", font="lato", palette=(
"#f15854", "#faa43a", "#e5d00d", "#60bd68", "#5da5da", "#f17cb0", "#975597", "#b2912f", "#aaaaaa", "#4d4d4d"))
# sns.set(palette=('#fcfbfd', '#efedf5', '#dadaeb', '#bcbddc', '#9e9ac8', '#807dba', '#6a51a3', '#54278f', '#3f007d'))
# sns.set(palette=(
# "#f15854", "#faa43a", "#e5d00d", "#60bd68", "#5da5da", "#f17cb0", "#975597", "#b2912f", "#aaaaaa", "#4d4d4d"))
# sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})
plt.rcParams['font.weight'] = "medium"
plt.rcParams['axes.labelweight'] = 'medium'
plt.rcParams['figure.titleweight'] = 'medium'
plt.rcParams['axes.titleweight'] = 'medium'
plt.rcParams['figure.figsize'] = 12, 3
plt.rcParams["legend.framealpha"] = 1
plt.rcParams['axes.xmargin'] = 0
plt.rcParams['axes.ymargin'] = 0

comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
get_color = {'A': '#0074FF', 'T': '#7900C8', 'G': '#FF4500', 'C': '#FFA500'}
nucs = ['A', 'C', 'G', 'T']
nuc_pos = dict(zip(nucs, range(4)))

field = 'fdrp_bh'
es_field = 'es_mean'

perf_tr = 0.0005
fdr_tr = 0.05

dpi = 300

text_h = -0.25*3/2
motive_pos_h = -0.5*4/3
text_color = '#505050'

# top10_names = ['CTCF_HUMAN.tsv', 'SPI1_HUMAN.tsv',
#                'FOXA1_HUMAN.tsv', 'CEBPB_HUMAN.tsv',
#                'ANDR_HUMAN.tsv', 'ESR1_HUMAN.tsv',
#                'NRF1_HUMAN.tsv', 'DUX4_HUMAN.tsv',
#                'CREB1_HUMAN.tsv', 'AP2A_HUMAN.tsv']
#
# motive_lengths = dict(zip(top10_names, [19, 17, 12, 12, 18, 15, 17, 11, 11, 15]))
# logos = dict(zip(top10_names, [name.replace('.tsv', '.H11MO.0.A.png') for name in top10_names]))

folder = '~/Releases/TF_P-values/TF_P-values/'
# folder = '~/DataForFigures/scatter/'
for name in os.listdir(os.path.expanduser(folder)):
    if 'Mathelier' in name or 'All' in name:
        continue
    if not os.path.isfile(os.path.expanduser("{}{}".format(folder, name))):
        continue
    print(name)
    pt = pd.read_table(os.path.expanduser("{}{}".format(folder, name)))
    try:
        pcm_file = os.path.expanduser("~/pcm/{}".format(name.replace('.tsv', '.H11MO.0.A.pcm')))
        with open(pcm_file, 'r') as f:
            motive_length = len(f.readlines()) - 1
    except FileNotFoundError:
        try:
            pcm_file = os.path.expanduser("~/pcm/{}".format(name.replace('.tsv', '.H11MO.0.B.pcm')))
            with open(pcm_file, 'r') as f:
                motive_length = len(f.readlines()) - 1
        except FileNotFoundError:
            try:
                pcm_file = os.path.expanduser("~/pcm/{}".format(name.replace('.tsv', '.H11MO.0.C.pcm')))
                with open(pcm_file, 'r') as f:
                    motive_length = len(f.readlines()) - 1
            except FileNotFoundError:
                print('Poop')
                continue
    logo_file = os.path.expanduser(pcm_file.replace('/pcm/', '/DataForFigures/logos/') + '.png')

    fsx = motive_length + 2.4
    fsy = 8

    left_pad = 1.4 / fsx
    right_pad = 1 / fsx
    top_pad = 0.05
    bottom_pad = 0.05
    image_h = 0.15
    image_bottom_pad = 0

    # print(pt.info())

    pt = pt[(pt['motif_log_pref'] >= -np.log10(perf_tr)) | (pt['motif_log_palt'] >= -np.log10(perf_tr))]
    pt = pt[~(pt[field + '_alt'].isnull() | pt[field + '_ref'].isnull())]
    pt = pt[(pt[field + '_alt'] <= fdr_tr) | (pt[field + '_ref'] <= fdr_tr)]

    pt['log_pv'] = (np.log10(
        pt[[field + '_ref', field + '_alt']]).min(axis=1)) \
                   * np.sign(pt[field + '_alt'] - pt[field + '_ref'])

    x = []
    y = []
    c = []
    for index, row in pt.iterrows():
        if row['log_pv'] < 0:
            major = row['ref']
            minor = row['alt']
            es = row[es_field + '_ref']
        elif row['log_pv'] > 0:
            major = row['alt']
            minor = row['ref']
            es = row[es_field + '_alt']

        if row['motif_orient'] == '-':
            major = comp[major]
            minor = comp[minor]

        x.append((row['motif_pos'] + (nuc_pos[minor] + 1) / 5) / motive_length)
        y.append(min(es / np.log(2), 3.98))  # FIXME
        c.append(get_color[major])

    if not x:
        print('empty')
        continue

    fig, ax = plt.subplots(dpi=dpi, figsize=(fsx, fsy))
    # fig.tight_layout(pad=0)
    ax.set_xticks([])
    ax.set_yticks([x/2 for x in range(0, 9)])

    ax.set_xlim(0, 1)
    ax.set_ylim(-0.8, 4)

    ax.set_position([left_pad, image_h + bottom_pad, 1 - left_pad - right_pad, 1 - image_h - top_pad - bottom_pad])

    im = Image.open(logo_file)
    im = im.resize((int(dpi * fsx * (1 - left_pad - right_pad)),
                    int(dpi * fsy * (image_h - image_bottom_pad))), Image.ANTIALIAS)
    IM = fig.figimage(im, dpi * fsx * left_pad, dpi * fsy * (image_bottom_pad + bottom_pad - 0.0125))

    plt.axhline(y=0, linewidth=2, linestyle='--', color='grey')
    for k in range(motive_length - 1):
        plt.axvline(x=(k + 1) / motive_length, linewidth=2, linestyle='--', color='grey')
    for k in range(motive_length):
        plt.text(x=(k + 0.5) / motive_length, y=motive_pos_h, s=str(k+1), color=text_color, ha='center', va='center')
        for i in range(4):
            plt.axvline(x=(k + (i + 1) / 5) / motive_length, linewidth=1, linestyle='--', color='grey')
            plt.text(x=(k + (i + 1) / 5) / motive_length, y=text_h, s=nucs[i], color=text_color,
                     ha='center', va='center', fontdict={'size': 18})
    plt.text(x=0.5 / motive_length, y=0.5 * text_h, s='minor', color=text_color,
             ha='center', va='center', fontdict={'size': 18})

    ax.scatter(x=x, y=y, c=c, s=30, zorder=10, alpha=0.7)
    plt.ylabel('Effect size, log₂')

    ax.text(x=0.99, y=3.95, s=name.replace('_HUMAN.tsv', ''), color=text_color, va='top', ha='right')

    plt.savefig(os.path.expanduser('~/AC_6/AS_Figure_6_stave_{}.svg'.format(name.replace('_HUMAN.tsv', ''))))
    plt.close(fig)
