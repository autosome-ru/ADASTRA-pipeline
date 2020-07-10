from svgutils import compose, transform
import numpy as np
from scipy.special import gamma
from math import lgamma
import sys
import os


def get_KDIC(counts, N):
    return 1 / (N * np.log(0.25)) * (
                lgamma(N + 1) + sum([x * np.log(0.25) - lgamma(x + 1) for x in counts]))


def get_heights(pcm_file, mode='freq'):
    heights = []
    lines = []
    with open(pcm_file) as f:
        for line in f:
            if line.startswith('>'):
                continue
            lines.append(list(map(float, line.strip('\n').split('\t'))))

    m = len(lines)
    for counts in lines:
        N = sum(counts)
        if mode == 'freq':
            heights.append(sorted(list(zip(letters, [x / N for x in counts])), key=lambda x: x[1], reverse=True))
        elif mode == 'KDIC':
            KDIC = get_KDIC(counts, N)
            print(KDIC)
            heights.append(sorted(list(zip(letters, [(x / N) * KDIC for x in counts])), key=lambda x: x[1], reverse=True))

    return m, heights


def place_letter_on_svg(figure, letter_svg, x, y, h):
    letter_object = transform.fromfile(letter_svg)
    letter_root = letter_object.getroot()
    # 13.229 and 26.458 are letter svg view box w and h
    letter_root.scale_xy(unit_width/13.229, h/26.458)
    letter_root.moveto(x, y)
    figure.append(letter_root)


def renorm(position):
    letters, heights = zip(*position)
    total_height = sum(heights)
    new_total_height = 0
    new_heights = []
    for height in heights:
        if height < visible_cut_tr * total_height:
            new_heights.append(0)
        else:
            new_total_height += height
            new_heights.append(height)
    new_heights = [x * new_total_height / total_height for x in new_heights]
    return zip(letters, new_heights)


if __name__ == '__main__':
    pcm_path = os.path.expanduser('~/pcm/CTCF_HUMAN.H11MO.0.A.pcm')

    letters = ['A', 'C', 'G', 'T']

    unit_width = 300
    unit_height = 600

    visible_cut_tr = 0.01

    revcomp = True

    letter_svgs = {
        'A': os.path.expanduser('~/letters/lettarA_path.svg'),
        'C': os.path.expanduser('~/letters/lettarC_path.svg'),
        'G': os.path.expanduser('~/letters/lettarG_path.svg'),
        'T': os.path.expanduser('~/letters/lettarT_path.svg'),
    }

    get_revcomp = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
    }

    m, heights = get_heights(pcm_path, mode='KDIC')
    fig = transform.SVGFigure("{}px".format(m * unit_width), "{}px".format(unit_height))

    for pos, pack in enumerate(heights[::-1] if revcomp else heights):
        current_height = 0
        for letter, height in renorm(pack):
            # Draw letter with offset of pos*unit_width, current_height*unit_height and height of height*unit_height
            place_letter_on_svg(fig, letter_svgs[get_revcomp[letter] if revcomp else letter], pos*unit_width, (1-current_height - height)*unit_height, height*unit_height)
            current_height += height

    fig.save(os.path.expanduser('~/test.svg'))
