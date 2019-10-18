import json
import math
import numpy as np
# import operator as op
# from functools import reduce
import os.path
import sys
# from multiprocessing import Process
# from matplotlib import pyplot as plt
# import random
import time

sys.path.insert(1, '/home/abramov/ASB-Project')

parameters_path = '/home/abramov/PARAMETERS/'
ploidy_path = "/home/abramov/PloidyForRelease/"


class Segmentation:
    def __init__(self):
        self.b = None
        self.bnum = None
        self.sc = None
        self.LINES = None
        self.length = None
        self.CRITICAL_GAP = None
        self.positions = None
        self.i_idxs = None
        self.last_snp_number = None
        self.candidate_numbers = None
        self.i_list = None
        self.candidates_count = None
        self.chrom = None
        self.mode = None
        self.dtype = np.float64
        self.bposn = None
        self.C = None
        self.S = None
        self.L = None
        self.P = None
        self.prior = None
        self.b_penalty = None
        self.d_penalty = None

    def loglikelyhood(self, N, X, i):
        p = 1 / (1 + i)
        if (self.mode == 'corrected' and N == 2 * X) or self.mode == 'binomial':
            return X * np.log(p) + (N - X) * np.log(1 - p) + np.log(self.prior[i])
        elif self.mode == 'corrected':
            return X * np.log(p) + (N - X) * np.log(1 - p) + np.log(self.prior[i]) + math.log(1 + i ** (2 * X - N))

    def get_P(self, first, last):
        if last - first == 1:
            return self.chrom.P_init[:, last]
        else:
            return np.sum(self.chrom.P_init[:, first + 1:last + 1], axis=1)

    def construct_probs(self):
        P = np.zeros((len(self.i_list), self.candidates_count + 1, self.candidates_count + 1), dtype=self.dtype)
        S = np.zeros((len(self.i_list), self.candidates_count + 1), dtype=self.dtype)
        self.L = np.zeros((self.candidates_count + 1, self.candidates_count + 1),
                          dtype=self.dtype)  # if in init -> -memory
        for j in range(0, self.candidates_count + 1):
            if j == 0:
                first = -1
            else:
                first = self.candidate_numbers[j - 1]
            if j == self.candidates_count:
                last = self.last_snp_number
            else:
                last = self.candidate_numbers[j]

            S[:, j] = self.get_P(first, last)
        self.S = S
        self.P = P

    def modify_P(self):
        # print('Constructing p-matrix')
        self.C = np.cumsum(self.S, axis=1)
        for j in range(self.candidates_count + 1):
            if j == 0:
                substract = 0
            else:
                substract = self.C[:, j - 1]
            for k in range(j, self.candidates_count + 1):
                self.P[:, j, k] = self.C[:, k] - substract

    def modify_L(self):
        # print('Constructing L')
        Q = np.sort(self.P, axis=0)
        self.L[:, :] = Q[-1, :, :] + np.log1p(np.sum(np.exp(Q[:-2, :, :] - Q[-1, :, :]), axis=0))

    def check_L(self):
        L = np.zeros((self.candidates_count + 1, self.candidates_count + 1), dtype=self.dtype)
        for j in range(self.candidates_count + 1):
            for k in range(j, self.candidates_count + 1):
                L[j, k] = 0
                accum = self.dtype(0)
                tmp = sorted([self.P[i, j, k] for i in self.i_idxs])
                for i in range(len(tmp) - 1):
                    accum += np.exp(tmp[i] - tmp[-1])
                L[j, k] += tmp[-1] + np.log1p(accum)
                if j != k and abs(-L[j, k] + self.L[j, k]) > 1:
                    print(j, k, L[j, k], round(-L[j, k] + self.L[j, k], 2))

    def get_distance_penalty(self, k):
        if self.d_penalty == 'NOTR' or \
                    (self.d_penalty == 'TR' and self.positions[k + 1] - self.positions[k] > self.CRITICAL_GAP):
            dist = (self.positions[k + 1] - self.positions[k]) / self.length
            return -1 * self.chrom.LINES * self.chrom.D_FACTOR * np.log1p(-dist)
        elif self.d_penalty == 'NONE' or self.d_penalty == 'TR':
            return 0
        else:
            raise ValueError(self.d_penalty)

    def get_parameter_penalty(self, borders, alphabet):
        k = borders * alphabet
        if self.b_penalty == 'CAIC':
            return -1 / 2 * k * (np.log(self.LINES) + 1)
        elif self.b_penalty == 'AIC':
            return -1 / 2 * k
        elif self.b_penalty == 'SQRT':
            return -1 / 2 * k * (np.sqrt(self.LINES) + 1)
        else:
            raise ValueError(self.b_penalty)
        
    def find_optimal_borders(self):
        # print('Constructing borders')
        for i in range(self.candidates_count + 1):
            self.sc[i] = self.L[0, i]

            kf = -1
            current_optimal = self.sc[i]

            for k in range(i):
                distance_penalty = self.get_distance_penalty(self.candidate_numbers[k] - self.candidate_numbers[0])
                parameter_penalty = self.get_parameter_penalty(self.bnum[k] + 1, len(self.i_list))

                likelyhood = self.sc[k] + self.L[k + 1, i] + distance_penalty
                candidate = likelyhood + parameter_penalty
                if candidate > current_optimal:
                    current_optimal = candidate
                    self.sc[i] = likelyhood
                    kf = k
            if kf != -1:
                self.b[kf] = True
            for j in range(kf + 1, i):
                self.b[j] = False

            self.bnum[i] = self.b.count(True)

        self.bposn = [self.candidate_numbers[i] for i in range(self.candidates_count) if self.b[i]]

    def estimate(self):
        self.construct_probs()
        self.modify_P()
        self.modify_L()
        # self.check_L()
        self.find_optimal_borders()


class PieceSegmentation(Segmentation):
    def __init__(self, chrom, start, end):
        super().__init__()
        self.start = start
        self.end = end
        self.CHR = chrom.CHR
        self.LINES = end - start + 1
        self.i_list = chrom.i_list
        self.i_idxs = range(len(self.i_list))
        self.positions = chrom.positions[
                         chrom.candidate_numbers[start]:chrom.candidate_numbers[end - 1] + 2]  # positions of snps
        self.candidate_numbers = chrom.candidate_numbers[start:end + 1]
        self.candidates_count = end - start
        if self.end == chrom.candidates_count:
            self.last_snp_number = chrom.LINES - 1
        else:
            self.last_snp_number = chrom.candidate_numbers[end + 1] - 1
        self.length = chrom.length  # (end - start + 1)/chrom.LINES*chrom.length
        self.CRITICAL_GAP = chrom.CRITICAL_GAP
        self.bposn = []  # border positions, snip numbers, border after
        self.sc = [0] * (self.candidates_count + 1)  # sc[i] = best log-likelyhood among all segmentations of snps[0,i]
        self.b = [False] * self.candidates_count  # bool borders, len=LINES.
        self.bnum = [0] * (self.candidates_count + 1)  # bnum[i] = number of borders before ith snp in best segmentation
        self.S = None
        self.C = None
        self.P = None  # segment-wise log-likelyhoods for each ploidy
        self.L = None
        self.chrom = chrom

        self.mode = chrom.mode  # binomial or corrected
        
        self.prior = chrom.prior
        self.b_penalty = chrom.b_penalty
        self.d_penalty = chrom.d_penalty

    # print(self.start, self.end, self.candidates_count, len(self.positions))


class ChromosomeSegmentation(Segmentation):  # chrom
    def __init__(self, seg, CHR, prior, i_list=None, length=0, n_max=0, CGF=0.0, d_penalty='CAIC', b_penalty='TR'):
        super().__init__()
        self.CHR = CHR  # name
        self.COV_TR = seg.COV_TR  # coverage treshold
        self.SEG_LENGTH = seg.SEG_LENGTH  # length of segment
        self.D_FACTOR = seg.D_FACTOR  # Multiplier in distance penalty
        self.INTERSECT = seg.INTERSECT  # length of intersection
        self.FILE = open(seg.FILE, 'r')
        self.SNPS, self.LINES = self.read_file_len()  # number of snps
        self.length = length  # length, bp
        if self.LINES == 0:
            return
        self.start = 0
        self.end = self.LINES - 2
        self.candidate_numbers = [i for i in range(self.LINES - 1)]
        self.candidates_count = self.LINES - 1
        self.last_snp_number = self.LINES - 1
        self.positions = []  # positions of snps
        self.ploidys = []  # max(ref, alt)/min(ref, alt)
        self.P = None  # segment-wise log-likelyhoods for each ploidy
        self.P_init = None  # snp-wise log0likelyhoods for each ploidy
        self.L = None  # segment-wise marginal log-likelyhoods
        self.CRITICAL_GAP = self.length * (
                1 - 10 ** (-CGF / self.LINES))  # critical distance after which border likelyhood increases
        self.bpos = []  # border positions between snps at x1 and x2: (x1+x2)/2 for x2-x1<=CRITICAL_GAP, (x1, x2) else
        self.bposn = []
        self.border_numbers = None
        self.chrom = self

        self.mode = seg.mode  # binomial or corrected

        if i_list:
            self.i_list = i_list
        else:
            self.i_list = []  # list of ploidys
            for i in range(1, n_max + 1):
                self.i_list.append(i)
        self.i_idxs = range(len(self.i_list))

        self.LS = None  # likelyhoods of splited segments for each ploidy
        self.ests = []  # estimated ploidys for splited segments
        self.quals = []  # qualities of estimations
        self.Q1 = []
        self.counts = []  # number of snps in segments
        self.second_ests = []
        
        self.prior = prior
        self.b_penalty = b_penalty
        self.d_penalty = d_penalty

    # [1017, 1035, 1160, 1442, 1529, 1641, 1857, 2045, 4062][2, 1, 5, 1, 3, 4, 2, 2, 1, 1]
    # self.set_candidates(set(range(5047)) - {1283, 1284})#set(random.sample(range(5047), 500)))
    # self.set_candidates({1017, 1035, 1160, 1283, 1442, 1641, 1857, 2045, 4062})# - {1283})
    # #| set([54, 114, 368, 979, 992, 1017, 1035, 1095, 1160, 1166, 1282, 1389, 1442, 1641,
    # 1725, 1842, 1857, 2045, 2071, 2710, 2797, 2798, 3030, 3223, 4062, 4345, 4362, 4586, 4589]) | {1283})#
    # | set(random.sample(range(5000), 20)))

    @staticmethod
    def get_params(line):
        line = line.split()
        chr = line[0]
        pos = int(line[1])
        ID = line[2]
        ref = line[3]
        alt = line[4]
        ref_c = int(line[5])
        alt_c = int(line[6])
        return chr, pos, ID, ref, alt, ref_c, alt_c

    def line_is_ok(self, line):
        if line[0] == '#':
            return False
        chr, pos, ID, ref, alt, ref_c, alt_c = self.get_params(line)
        if chr != self.CHR or ID == '.':
            return False
        return pos, ref_c, alt_c

    def read_file_len(self):
        count = 0
        snps = []
        for line in self.FILE:
            if not self.line_is_ok(line):
                continue
            count += 1
            snps.append(line)
        self.FILE.seek(0)
        return snps, count

    @staticmethod
    def split_list(length, l, k):
        iterator = []
        if length < l:
            iterator.append((0, length - 1))
            return iterator
        length -= k
        div, mod = divmod(length - 1, l - k)
        new_l, num = divmod(length - 1, div)
        for i in range(div):
            if i < num:
                iterator.append(((new_l + 1) * i, (new_l + 1) * (i + 1) + k))
            else:
                iterator.append((new_l * i + num, new_l * (i + 1) + k + num))
        return iterator

    def set_candidates(self, candidate_set):
        self.candidate_numbers = sorted(list(candidate_set))
        self.candidates_count = len(self.candidate_numbers)
        self.end = self.candidates_count - 1

    def construct_probs_initial(self):
        current_snip = -1

        S = np.zeros((len(self.i_list), self.LINES), dtype=self.dtype)
        for j in range(0, self.LINES):
            line = self.chrom.SNPS[self.start + j]

            result = self.chrom.line_is_ok(line)

            if not result:
                continue

            pos, ref_c, alt_c = result
            current_snip += 1
            N = ref_c + alt_c
            X = min(ref_c, alt_c)

            self.positions.append(pos)

            if X == 0:
                self.ploidys.append(0)
            else:
                self.ploidys.append(N / X - 1)

            for i in self.i_idxs:
                assert (self.i_list[i] > 0)
                S[i, current_snip] = self.loglikelyhood(N, X, self.i_list[i])
        self.P_init = S

    def find_optimal_borders(self):
        self.sc = [0] * (self.candidates_count + 1)  # sc[i] = best log-likelyhood among all segmentations of snps[0,i]
        self.b = [(0, 0)] * self.candidates_count  # borders, len=LINES. b[i]: (0,0) if there is no border after ith snp
        # (1, x) else,
        # where x is 1 if positions[i+1]-positions[i] > CRITICAL_GAP, 0 else
        self.bnum = [0] * (self.candidates_count + 1)  # bnum[i] = number of borders before ith snp in best segmentation
        for i in range(self.candidates_count + 1):
            self.sc[i] = self.L[0, i]

            kf = -1
            current_optimal = self.sc[i]

            for k in range(i):

                distance_penalty = self.get_distance_penalty(self.candidate_numbers[k])
                parameter_penalty = self.get_parameter_penalty(self.bnum[k] + 1, len(self.i_list))

                likelyhood = self.sc[k] + self.L[k + 1, i] + distance_penalty
                candidate = likelyhood + parameter_penalty
                if candidate > current_optimal:
                    current_optimal = candidate
                    self.sc[i] = likelyhood
                    kf = k
            if kf != -1:
                if self.positions[self.candidate_numbers[kf] + 1] - self.positions[self.candidate_numbers[kf]] \
                        > self.CRITICAL_GAP:
                    self.b[kf] = (1, 1)
                else:
                    self.b[kf] = (1, 0)
            for j in range(kf + 1, i):
                self.b[j] = (0, 0)

            self.bnum[i] = [x[0] for x in self.b].count(1)
            assert ([x[0] for x in self.b].count(1) == self.bnum[i])

        self.bposn = [self.candidate_numbers[i] for i in range(self.candidates_count) if self.b[i][0]]
        self.border_numbers = [-1] + [i for i in range(self.candidates_count) if self.b[i][0]] + [self.candidates_count]
        acum_counts = [0] + self.bposn + [self.last_snp_number]
        self.counts = [acum_counts[i + 1] - acum_counts[i] for i in range(len(acum_counts) - 1)]

        for i in range(len(self.b)):
            if self.b[i][0]:
                if self.b[i][1]:
                    self.bpos.append(
                        (self.positions[self.candidate_numbers[i]], self.positions[self.candidate_numbers[i] + 1]))
                else:
                    self.bpos.append(
                        (self.positions[self.candidate_numbers[i]] + self.positions[self.candidate_numbers[i] + 1]) / 2)

    def estimate_Is(self):
        self.LS = np.zeros(len(self.i_list), dtype=self.dtype)
        for n in range(len(self.border_numbers) - 1):
            first = self.border_numbers[n] + 1
            last = self.border_numbers[n + 1]
            self.LS = self.P[:, first, last] - self.L[first, last]
            i_max = np.argmax(self.LS)
            if i_max == 0:
                left_qual = 0
                Q1 = -1 * int(round(self.LS[1] - self.LS[0]))
            else:
                left_qual = -1 * int(round(self.LS[i_max - 1] - self.LS[i_max]))
                Q1 = -1 * int(round(self.LS[i_max] - self.LS[0]))
            if i_max == len(self.i_list) - 1:
                right_qual = 0
            else:
                right_qual = -1 * int(round(self.LS[i_max + 1] - self.LS[i_max]))
            # noinspection PyTypeChecker
            self.ests.append(self.i_list[i_max])
            self.quals.append((left_qual, right_qual))
            self.Q1.append(Q1)

        print('Total SNPs: {}, estimated ploidys: {}, border positions: {}'.format(len(self.positions), self.ests,
                                                                                   self.bpos))

    def estimate_chr(self):
        start_t = time.clock()

        if self.LINES == 0:
            return

        self.construct_probs_initial()

        if self.candidates_count > self.SEG_LENGTH:
            tuples = self.split_list(self.candidates_count + 1, self.SEG_LENGTH, self.INTERSECT)
            border_set = set()
            print("{} segments: {}".format(len(tuples), tuples))
            counter = 0
            for first, last in tuples:
                counter += 1
                print('Making {} out of {} segments from {} to {} for {}.'.format(counter, len(tuples), first, last,
                                                                                  self.CHR))
                PS = PieceSegmentation(self, first, last)
                PS.estimate()
                # print(PS.bposn)
                border_set |= set(PS.bposn)
            self.candidate_numbers = sorted(list(border_set))
            self.candidates_count = len(self.candidate_numbers)
        print(len(self.positions))
        print('{} candidates'.format(self.candidates_count))
        #		print('{} candidates: {}'.format(self.candidates_count, [x for x in self.candidate_numbers]))

        self.estimate()
        self.estimate_Is()
        print(self.bposn, self.ests, self.counts)
        print('CHR time: {} s'.format(time.clock() - start_t))

        # self.plot()

    def score_segmentation(self, bposn, ests):
        if len(ests) - len(bposn) != 1:
            print('Invalid segmentation')
            return 0
        score = 0
        bposn = [-1] + bposn + [self.LINES - 1]
        for i in range(len(bposn) - 1):
            first = bposn[i] + 1
            last = bposn[i + 1] + 1
            score += sum([self.P_init[ests[i]][j] for j in range(first, last)])
            for k in bposn[1:-1]:
                score += self.get_distance_penalty(k)
            score += self.get_parameter_penalty(len(bposn) - 2, len(self.i_list))
        return score


class GenomeSegmentator:  # seg
    def __init__(self, file, out, segm_mode, extra_states=None, b_penalty='CAIC', d_penalty='TR'):
        chr_l = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
                 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
                 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468,
                 156040895, 57227415]
        self.chrs, self.chr_lengths = self.construct_chrs(chr_l, from_sys=False)

        self.mode = segm_mode
        if extra_states:
            self.i_list = [1, 2, 3, 4, 5] + extra_states
        else:
            self.i_list = [1, 2, 3, 4, 5]

        self.FILE = file  # table
        self.OUT = open(out, 'w')  # ploidy file
        self.n_max = 5  # max ploidy
        self.D_FACTOR = 0.25
        self.CRITICAL_GAP_FACTOR = 16.5
        self.NUM_TR = 100  # minimal number of snps in chromosome to start segmentation
        self.COV_TR = 0  # coverage treshold
        self.INTERSECT = 300
        self.SEG_LENGTH = 600
        self.ISOLATED_SNP_FILTER = 4
        self.chr_segmentations = []  # chroms
        
        self.b_penalty = b_penalty
        self.d_penalty = d_penalty
        self.prior = dict(zip(self.i_list, [1]*len(self.i_list)))

        for CHR in self.chrs:
            chrom = ChromosomeSegmentation(self, CHR, self.prior, self.i_list, self.chr_lengths[CHR], self.n_max,
                                           self.CRITICAL_GAP_FACTOR, b_penalty=self.b_penalty, d_penalty=self.d_penalty)
            print('{} total SNP count: {}'.format(CHR, chrom.LINES))
            self.chr_segmentations.append(chrom)

    @staticmethod
    def construct_chrs(chr_l, from_sys=False):
        chrs = []
        chr_lengths = dict()
        if from_sys:
            chrs = sys.argv[3:]
        else:
            for i in range(1, 23):
                chrs.append('chr' + str(i))
            chrs.append('chrX')
            chrs.append('chrY')

        for i in range(1, 23):
            chr_lengths['chr' + str(i)] = chr_l[i - 1]
        chr_lengths['chrX'] = chr_l[-2]
        chr_lengths['chrY'] = chr_l[-1]

        return chrs, chr_lengths

    def append_ploidy_segments(self, chrom):
        segments_to_write = []
        cur = 1
        counter = 0
        if chrom.LINES >= self.NUM_TR:
            for border in chrom.bpos + [chrom.length]:
                if isinstance(border, tuple):
                    segments_to_write.append([chrom.CHR, cur, border[0] + 1, chrom.ests[counter], chrom.Q1[counter],
                                              chrom.quals[counter][0], chrom.quals[counter][1],
                                              chrom.counts[counter]])
                    cur = border[0] + 1
                    segments_to_write.append([chrom.CHR, cur, border[1], 0, 0, 0, 0, 0])
                    cur = border[1]
                    counter += 1
                else:
                    segments_to_write.append(
                        [chrom.CHR, cur, math.floor(border) + 1, chrom.ests[counter], chrom.Q1[counter],
                         chrom.quals[counter][0], chrom.quals[counter][1],
                         chrom.counts[counter]])
                    cur = math.floor(border) + 1
                    counter += 1
        else:
            segments_to_write.append([chrom.CHR, 1, chrom.length, 0, 0, 0, 0, 0])
        return segments_to_write

    def write_ploidy_to_file(self, chrom):
        segments = self.append_ploidy_segments(chrom)

        filtered_segments = self.filter_segments(segments, self.ISOLATED_SNP_FILTER)
        for segment in filtered_segments:
            self.OUT.write('\t'.join(map(str, segment)) + '\n')

    # noinspection PyTypeChecker
    def estimate_ploidy(self):
        for j in range(len(self.chr_segmentations)):
            chrom = self.chr_segmentations[j]
            chrom.estimate_chr()
            self.write_ploidy_to_file(chrom)
            self.chr_segmentations[j] = None

    @staticmethod
    def filter_segments(segments, snp_number_tr=3):
        is_bad_left = False
        bad_segments_indexes = set()
        is_bad_segment = False
        for k in range(len(segments)):
            if segments[k][7] <= snp_number_tr:  # если k сегмент "плохой"
                if is_bad_segment:  # если k-1 тоже "плохой"
                    bad_segments_indexes.add(k - 1)  # убрать k-1
                    is_bad_left = True
                else:
                    is_bad_left = False
                is_bad_segment = True  # текущий сегмент плохой, следующий шаг цикла
            else:  # k сегмент хороший
                if is_bad_segment and not is_bad_left and k > 1:  # а k-1 плохой и k-2 хороший
                    if segments[k][3] < segments[k - 1][3] and segments[k - 2][3] < segments[k - 1][3]:
                        # если CNR k-1 сегмента больше CNR k-2 и k сегментов
                        if segments[k][3] > segments[k - 2][3]:  # если CNR k сегмента больше CNR k-2
                            segments[k - 1][3] = segments[k][3]  # присвоить CNR k сегмента
                        elif segments[k][3] < segments[k - 2][3]:  # если CNR k-2 сегмента больше CNR k
                            segments[k - 1][3] = segments[k - 2][3]  # присвоить CNR k-2 сегмента
                        else:  # CNR равны
                            segments[k - 1][3] = segments[k][3]

                        for j in range(4, 7):
                            segments[k-1][j] = 0
                    is_bad_left = True
                is_bad_segment = False  # текущий сегмент хороший, следующий шаг цикла

        filtered_segments = []
        for k in range(len(segments)):
            if k in bad_segments_indexes:
                continue
            filtered_segments.append(segments[k])
        return filtered_segments


if __name__ == '__main__':
    key = sys.argv[1]
    print(key)
    
    mode = sys.argv[2].lower()
    states = [] if sys.argv[3] == '_' else list(map(float, sys.argv[3].split('_')))
    b_penalty = sys.argv[4]
    d_penalty = sys.argv[5]
    
    model = '-'.join(sys.argv[2:6])
    print(model)
    
    out_file = ploidy_path + key + ".tsv"

    t = time.clock()
    if not os.path.isdir(ploidy_path + model):
        os.mkdir(ploidy_path + model)
    GS = GenomeSegmentator(out_file, ploidy_path + model + key + "_ploidy.tsv", mode, states, b_penalty, d_penalty)
    GS.estimate_ploidy()
    print('Total time: {} s'.format(time.clock() - t))
