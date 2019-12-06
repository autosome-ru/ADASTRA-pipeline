import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from scipy import optimize
import operator as op
from functools import reduce


def ncr(n, r):
    r = min(r, n - r)
    numer = reduce(op.mul, range(n, n - r, -1), 1)
    denom = reduce(op.mul, range(1, r + 1), 1)
    return numer / denom


def make_ncr_array(n_max):
    rv = np.zeros((n_max + 1, n_max + 1), dtype=np.float128)
    for n in range(n_max + 1):
        n_pow = 2**(-n)
        for k in range(n + 1):
            rv[n, k] = ncr(n, k) * n_pow
    return rv


def make_counts_array(stats, n_max):
    c = np.zeros((n_max + 1, n_max + 1), dtype=np.int64)
    for n in range(n_max + 1):
        for k in range(n + 1):
            slice = stats[(stats['cover'] == n) & (stats['ref_counts'] == k)]
            if not slice.empty:
                assert len(slice.index) == 1
                c[n, k] = np.int64(slice['ref_counts'])
    return c


def make_target(counts, nck, n_min, n_max):
    def target(alpha):
        return -1 * sum(
            sum(counts[n, k] * np.log(nck[n, k] * (1 - alpha) + 2 * k / (n * (n + 1)) * alpha)
                for k in range(0, n + 1))
            for n in range(n_min, n_max + 1))

    return target


def make_derivative(counts, nck, n_min, n_max):
    def target(alpha):
        return -1 * sum(
            sum(counts[n, k] * (nck[n, k] * (-1) + 2 * k / (n * (n + 1))) /
                (nck[n, k] * (1 - alpha) + 2 * k / (n * (n + 1)) * alpha)
                for k in range(0, n + 1))
            for n in range(n_min, n_max + 1))

    return target


if __name__ == '__main__':
    n_min = 10
    n_max = 50
    stats = pd.read_table('~/cover_bias_statistics.tsv')
    counts = make_counts_array(stats, n_max)
    print('made counts')
    nck = make_ncr_array(n_max)
    print('made ncr')

    x = [(a/1000) for a in range(1, 101)]
    values = [make_derivative(counts, nck, n_min, n_max)(a/100) for a in range(1, 101)]

    plt.scatter(x, values)
    plt.grid(True)


    print(optimize.brenth(f=make_derivative(counts, nck, n_min, n_max), a=0, b=0.99))

    # print(optimize.minimize(fun=make_target(counts, nck, n_min, n_max), x0=0, method='TNC',
    #                         jac=make_derivative(counts, nck, n_min, n_max), bounds=[(0, 0.5)]))
    plt.show()