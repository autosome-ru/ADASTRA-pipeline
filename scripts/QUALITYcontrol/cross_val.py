import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns


class ExpectationCounter:
    def __init__(self, tr):
        self.tr = tr
        self.dict = {}

    def expected(self, n, p):
        if (n, p) in self.dict:
            return self.dict[(n, p)]
        else:
            return self._expected(n ,p)

    def _expected(self, n, p):
        dist = stats.binom(n=n, p=p)
        norm = dist.cdf(n-self.tr) - dist.cdf(self.tr-1)
        pmf = lambda x: dist.pmf(x) / norm
        value = sum(pmf(x)*x for x in range(self.tr, n//2)) + sum(pmf(x)*(n-x) for x in range(n//2, n - self.tr + 1))
        self.dict[(n, p)] = value
        return value


def chisq(arg):
    n, x, p = arg
    return (x-n*p)**2/(n*p*(1-p))


def g(arg, es):
    n, x, p = arg
    exp = es.expected(n, p)
    print(exp)
    return 2 * x * np.log(x / exp)


def get_p_value(arg):
    n, x, p = arg
    dist = stats.binom(n=n, p=p)
    cdf = dist.cdf
    return (cdf(x) - cdf(4) + cdf(n-5) - cdf(n-x-1)) / (cdf(n-5) - cdf(4)) if x < n/2 else 1


def get_lik(arg):
    n, x, p = arg
    dist = stats.binom(n=n, p=p)
    pmf = dist.pmf
    cdf = dist.cdf
    return (pmf(x) + (pmf(n-x) if x < n/2 else 0)) / (cdf(n-5) - cdf(4))


def calculate_gof(nxp_array, ec):
    observed = np.array([x for n, x, p in nxp_array], dtype=np.int64)
    expected = np.array([ec.expected(n, p) for n, x, p in nxp_array], dtype=np.float_)

    df = len(observed) - 1
    g_stat = np.sum(observed * np.log(observed / expected)) * 2
    norm = sum(observed)

    score = np.sqrt(max(g_stat - df, 0) / (df * (norm - 1)))

    return score


if __name__ == '__main__':
    N = 1000
    iters = 10

    summ = 0

    ec = ExpectationCounter(5)

    for iter in range(iters):
        cov = 10
        n = np.random.negative_binomial(cov*2, 0.5, N) + 1
        BAD = 1
        p = np.array([1/(BAD+1)]*N)
        BAD_test = 2
        p_test = np.array([1/(BAD_test+1)]*N)
        x = np.random.binomial(n, p)
        x = np.minimum(x, n-x)

        nxp = [(n, x, p) for n, x, p in zip(n, x, p_test) if 5 <= x <= n - 5]

        N = len(nxp)

        st = np.array(list(map(lambda f: g(f, ec), nxp)))
        s = np.sum(st)
        print(s, -np.log10(stats.distributions.chi2.sf(s, N)))
        summ += s

        pvs = np.array(list(map(get_p_value, nxp)))
        plt.hist(pvs, range=(0, 1), bins=50)


        rmsea = calculate_gof(nxp, ec)
        print(rmsea)
        plt.title('Sim: {:.2f}, Test: {:.2f}, RMSEA: {:.3f}\nCov: {}, SNPs: {}'.format(BAD, BAD_test, rmsea, cov, N))

        plt.show()

    print(summ/iters)
