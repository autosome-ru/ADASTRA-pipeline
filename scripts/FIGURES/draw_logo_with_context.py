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
            heights.append(sorted(list(zip(letters, [(x / N) * KDIC for x in counts])), key=lambda x: x[1], reverse=True))

    return m, heights


def place_drawing_on_svg(figure, drawing_svg, x, y, h, w, self_w, self_h):
    drawing_object = transform.fromfile(drawing_svg)
    drawing_root = drawing_object.getroot()
    drawing_root.scale_xy(w/self_w, h/self_h)
    drawing_root.moveto(x, y)
    figure.append(drawing_root)


def place_letter_on_svg(figure, letter_svg, x, y, h, w):
    # 13.229 and 26.458 are letter svg view box w and h
    place_drawing_on_svg(figure, letter_svg, x, y, h, w, 13.229, 26.458)


def place_dna_on_svg(figure, dna_svg, x, y, h=100, w=365):
    place_drawing_on_svg(figure, dna_svg, x, y, h, w, 365, 100)


def place_hill_on_svg(figure, hill_svg, x, y, h, w):
    place_drawing_on_svg(figure, hill_svg, x, y, h, w, 2420.548, 479.477)


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
    context = ' ' * 20 + 'gggcaacagggctcgggcggccccCggtggtcagggcctgcatgttggt' + ' ' * 20
    alt = 'T'
    pos_in_motif = 8
    ef = 1.13
    asb_is_ref = False
    revcomp = True

    letters = ['A', 'C', 'G', 'T']

    unit_width = 300
    unit_height = 600

    motif_gap = 0.1
    snp_gap = 0.1
    text_h = 0.5
    text_width = 0.85
    snp_text_h = 0.8
    snp_text_width = 1
    hill_sum_height = (2*snp_text_h + snp_gap - text_h)*1.1
    hill_width = 7
    pseudocount = 0.1
    add_letters = 6

    visible_cut_tr = 0.05

    indent = max(snp_text_h + (snp_gap - text_width) / 2, hill_sum_height * (1 - pseudocount))

    full_gap = motif_gap + indent

    letter_svgs = {
        'A': os.path.expanduser('~/letters/lettarA_path.svg'),
        'C': os.path.expanduser('~/letters/lettarC_path.svg'),
        'G': os.path.expanduser('~/letters/lettarG_path.svg'),
        'T': os.path.expanduser('~/letters/lettarT_path.svg'),
    }

    black_letter_svgs = {
        'A': os.path.expanduser('~/letters/lettarA_path_black.svg'),
        'C': os.path.expanduser('~/letters/lettarC_path_black.svg'),
        'G': os.path.expanduser('~/letters/lettarG_path_black.svg'),
        'T': os.path.expanduser('~/letters/lettarT_path_black.svg'),
    }

    hill_svgs = {
        'A': os.path.expanduser('~/letters/hill_A.svg'),
        'C': os.path.expanduser('~/letters/hill_C.svg'),
        'G': os.path.expanduser('~/letters/hill_G.svg'),
        'T': os.path.expanduser('~/letters/hill_T.svg'),
    }

    get_revcomp = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
    }

    m, heights = get_heights(pcm_path, mode='KDIC')
    print(full_gap, text_h, indent)
    fig = transform.SVGFigure("{}px".format((m + add_letters + add_letters) * unit_width), "{}px".format(unit_height * (1 + full_gap + text_h + indent)))

    pos_in_motif = m - pos_in_motif - 1 if revcomp else pos_in_motif
    motif_context = context[44 - add_letters - pos_in_motif: 44 - pos_in_motif + m + add_letters]
    pos_in_motif += add_letters

    hill_height = max(min(ef, 2) / 2, pseudocount)
    hill_height = (hill_height + 1) / 2
    ref_height = hill_sum_height * hill_height
    alt_height = hill_sum_height * (1 - hill_height)

    print(motif_context, pos_in_motif)

    place_letter_on_svg(fig, os.path.expanduser('~/letters/rect.svg'), pos_in_motif * unit_width, 0, (1 + full_gap + (alt_height - ref_height + text_h)/2 + snp_gap/2 + 2*snp_text_h) * unit_height, unit_width)

    for pos, pack in enumerate(heights[::-1] if revcomp else heights):
        pos += add_letters
        current_height = 0
        for letter, height in renorm(pack):
            # Draw letter with offset of pos*unit_width, current_height*unit_height and height of height*unit_height
            place_letter_on_svg(fig, letter_svgs[get_revcomp[letter] if revcomp else letter], pos*unit_width, (1-current_height - height)*unit_height, height*unit_height, unit_width)
            current_height += height

    for pos, letter in enumerate(motif_context):
        if letter == ' ':
            continue
        if pos != pos_in_motif:
            place_letter_on_svg(fig, black_letter_svgs[letter.upper()], (pos + (1 - text_width) / 2) * unit_width, (1 + full_gap) * unit_height, text_h * unit_height, text_width * unit_width)
        else:
            if not asb_is_ref:
                ref_height, alt_height = alt_height, ref_height
            place_hill_on_svg(fig, hill_svgs[letter], (pos - (hill_width - 1) / 2) * unit_width, (1 + full_gap - ref_height) * unit_height, ref_height * unit_height, hill_width * unit_width)
            place_hill_on_svg(fig, hill_svgs[alt], (pos - (hill_width - 1) / 2) * unit_width, (1 + full_gap + text_h + alt_height) * unit_height, -alt_height * unit_height, hill_width * unit_width)
            place_letter_on_svg(fig, letter_svgs[letter], (pos + (1 - snp_text_width) / 2) * unit_width, (1 + full_gap + (alt_height - ref_height + text_h)/2 - snp_text_h - snp_gap/2) * unit_height, snp_text_h * unit_height, snp_text_width * unit_width)
            place_letter_on_svg(fig, letter_svgs[alt], (pos + (1 - snp_text_width) / 2) * unit_width, (1 + full_gap + (alt_height - ref_height + text_h)/2 + snp_gap/2) * unit_height, snp_text_h * unit_height, snp_text_width * unit_width)

    text_x = pos_in_motif + 2
    txt_ref = transform.TextElement(text_x * unit_width, (1 + motif_gap + indent/2) * unit_height, 'P-value: 123', size=200)
    txt_alt = transform.TextElement(text_x * unit_width, (1 + motif_gap + text_h + indent*3/2) * unit_height, 'P-value: 123', size=200)
    fig.append([txt_ref, txt_alt])

    fig.save(os.path.expanduser('~/test.svg'))
