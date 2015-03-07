from __future__ import division
from macroeco_distributions import pln, pln_solver, negbin_solver
from numpy import empty
from scipy.stats import nbinom

"""Obtain ppf of a distribution function"""

def get_rad_pln(S, mu, sigma, lower_trunc = True):
    """Obtain the predicted RAD from a Poisson lognormal distribution"""
    abundance = list(empty([S]))
    rank = range(1, int(S) + 1)
    cdf_obs = [(rank[i]-0.5) / S for i in range(0, int(S))]
    j = 0
    cdf_cum = 0
    i = 1
    while j < S:
        cdf_cum += pln.pmf(i, mu, sigma, lower_trunc)
        while cdf_cum >= cdf_obs[j]:
            abundance[j] = i
            j += 1
            if j == S:
                abundance.reverse()
                return abundance
        i += 1
        
def get_rad_negbin(S, n, p):
    """Obtain the predicted RAD from a negative binomial distribution"""
    abundance = list(empty([S]))
    rank = range(1, int(S) + 1)
    cdf_obs = [(rank[i]-0.5) / S for i in range(0, int(S))]
    j = 0
    cdf_cum = 0
    i = 1
    while j < S:
        cdf_cum += nbinom.pmf(i, n, p) / (1 - nbinom.pmf(0, n, p))
        while cdf_cum >= cdf_obs[j]:
            abundance[j] = i
            j += 1
            if j == S:
                abundance.reverse()
                return abundance
        i += 1

def get_rad_from_obs(ab, dist):
    if dist == 'negbin':
        n, p = negbin_solver(ab)
        pred_rad = get_rad_negbin(len(ab), n, p)
    elif dist == 'pln':
        mu, sigma = pln_solver(ab)
        pred_rad = get_rad_pln(len(ab), mu, sigma)
    return pred_rad