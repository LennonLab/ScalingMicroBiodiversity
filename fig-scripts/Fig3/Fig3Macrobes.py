from __future__ import division
import  matplotlib.pyplot as plt
import numpy as np
import random
import scipy as sc
from scipy import stats
import os
import sys
import statsmodels.stats.api as sms
import statsmodels.formula.api as smf
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from statsmodels.stats.outliers_influence import summary_table
from numpy import log, log2, exp, sqrt, log10, pi
from scipy.optimize import fsolve
import scipy.optimize as opt
import pandas as pd #import patsy
import mpmath as mpm
from scipy.optimize import fsolve
from math import erf, pi
import linecache
import math

mydir = os.path.expanduser("~/GitHub/rare-bio/")
mydir2 = os.path.expanduser("~/")


def alpha2(a, N, Nmax, Nmin):
    y = sqrt(pi*Nmin*Nmax)/(2.0*a) * exp((a * log2(sqrt(Nmax/Nmin)))**2.0)
    y = y * exp((log(2.0)/(2.0*a))**2.0)
    y = y * erf(a * log2(sqrt(Nmax/Nmin)) - log(2.0)/(2.0*a)) + erf(a * log2(sqrt(Nmax/Nmin)) + log(2.0)/(2.0*a))
    y -= N

    return y # find alpha

def s2(a, Nmax, Nmin):
    return sqrt(pi)/a * exp( (a * log2(sqrt(Nmax/Nmin)))**2) # Using equation 10


def getNmax(N, b, slope):
    return 10 ** (b + slope*(log10(N)))

def expS(N, b, slope):
    return 10 ** (b + slope*(log10(N)))



def Fig3():

    """ A figure demonstrating a strong richness relationship across 10 or 11
    orders of magnitude in total abundance. Taxonomic richness of a sample
    scales in a log-log fashion with the total abundance of the sample.
    """

    Nlist, Slist, klist, NmaxList, datasets, radDATA = [[],[],[],[],[],[]]

    name = 'MCDB'
    kind = 'macro'
    path = mydir2+'data/macro/'+name+'/'+name+'-SADMetricData.txt'
    numlines = sum(1 for line in open(path))
    lines = np.random.choice(range(1, numlines+1), 1000, replace=True)

    path = mydir2+'data/'+kind+'/'+name+'/'+name+'-SADMetricData.txt'
    for line in lines:
        data = linecache.getline(path, line)
        radDATA.append(data)

    for data in radDATA:
        data = data.split()
        name, kind, N, S, Var, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = data

        N = float(N)
        S = float(S)
        Nmax = float(Nmax)

        if S < 2 or N < 11: continue # Min species richness

        Nlist.append(float(np.log10(N)))
        Slist.append(float(np.log10(S)))

        NmaxList.append(float(np.log10(Nmax)))
        klist.append('DarkCyan')


    Nlist, Slist, NmaxList = zip(*sorted(zip(Nlist, Slist, NmaxList)))
    Nlist = list(Nlist)
    Slist = list(Slist)
    NmaxList = list(NmaxList)

    # Regression for Dominance (Nmax) vs. N
    d = pd.DataFrame({'N': Nlist})
    d['Nmax'] = NmaxList
    f = smf.ols('Nmax ~ N', d).fit()

    dR2 = f.rsquared
    dintercept = f.params[0]
    dslope = f.params[1]
    print 'R2 for Nmax vs. N:', round(dR2,3)

    # Regression for Richness (S) vs. N
    d = pd.DataFrame({'N': Nlist})
    d['S'] = Slist
    f = smf.ols('S ~ N', d).fit()

    sR2 = f.rsquared
    sintercept = f.params[0]
    sslope = f.params[1]
    print 'R2 for S vs. N:', round(sR2,3),'\n'

    print 'Nmax =', round(dintercept,2), '*', 'N^', round(dslope,2)
    print 'S =', round(sintercept,2), '*', 'N^', round(sslope,2),'\n'

    MammalN = 286 * (10**6)
    MammalNmax = 82 * (10**6)
    MammalS = 63

    # UK estimates for mammals
    Nmin = 1
    N = float(MammalN)
    empS = expS(N, sintercept, sslope)

    Nmax = MammalNmax
    guess = 0.183774388267

    a = opt.fsolve(alpha2, guess, (N, Nmax, Nmin))[0]
    print guess, a
    S2 = s2(a, Nmax, Nmin)

    print 'MammalNmax:', '%.2e' % float(MammalNmax)
    print 'predicted Nmax:', '%.2e' % getNmax(N, dintercept, dslope),'\n'
    print 'estimated S for Mammals:', '%.2e' % MammalS
    print 'scaling law prediction of S for Mammals:', '%.2e' % empS
    print 'lognormal prediction of S for Mammals:', '%.2e' % S2,'\n'



    AvianN = 2.*10**11
    AvianNmax = 3*10**9
    AvianS = 10500

    Nlist, Slist, klist, NmaxList, datasets, radDATA = [[],[],[],[],[],[]]

    name = 'BBS'
    kind = 'macro'
    path = mydir2+'data/macro/'+name+'/'+name+'-SADMetricData.txt'
    numlines = sum(1 for line in open(path))
    lines = np.random.choice(range(1, numlines+1), 1000, replace=True)

    path = mydir2+'data/'+kind+'/'+name+'/'+name+'-SADMetricData.txt'
    for line in lines:
        data = linecache.getline(path, line)
        radDATA.append(data)

    for data in radDATA:
        data = data.split()
        name, kind, N, S, Var, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = data

        N = float(N)
        S = float(S)
        Nmax = float(Nmax)

        if S < 2 or N < 11: continue # Min species richness

        Nlist.append(float(np.log10(N)))
        Slist.append(float(np.log10(S)))

        NmaxList.append(float(np.log10(Nmax)))
        klist.append('DarkCyan')


    Nlist, Slist, NmaxList = zip(*sorted(zip(Nlist, Slist, NmaxList)))
    Nlist = list(Nlist)
    Slist = list(Slist)
    NmaxList = list(NmaxList)

    Nlist = list(Nlist)
    Slist = list(Slist)
    NmaxList = list(NmaxList)

    # Regression for Dominance (Nmax) vs. N
    d = pd.DataFrame({'N': Nlist})
    d['Nmax'] = NmaxList
    f = smf.ols('Nmax ~ N', d).fit()

    dR2 = f.rsquared
    dintercept = f.params[0]
    dslope = f.params[1]
    print 'R2 for Nmax vs. N:', round(dR2,3)

    # Regression for Richness (S) vs. N
    d = pd.DataFrame({'N': Nlist})
    d['S'] = Slist
    f = smf.ols('S ~ N', d).fit()

    sR2 = f.rsquared
    sintercept = f.params[0]
    sslope = f.params[1]
    print 'R2 for S vs. N:', round(sR2,3),'\n'

    print 'Nmax =', round(dintercept,2), '*', 'N^', round(dslope,2)
    print 'S =', round(sintercept,2), '*', 'N^', round(sslope,2),'\n'

    # Global estimates for birds
    Nmin = 1
    N = float(AvianN)
    empS = expS(N, sintercept, sslope)

    Nmax = AvianNmax
    guess = 0.210521390628


    a = opt.fsolve(alpha2, guess, (N, Nmax, Nmin))[0]
    print guess, a
    S2 = s2(a, Nmax, Nmin)

    print 'AvianNmax:', '%.2e' % float(AvianNmax)
    print 'predicted Nmax:', '%.2e' % getNmax(N, dintercept, dslope),'\n'
    print 'estimated S for Birds:', '%.2e' % AvianS
    print 'scaling law prediction of S for Birds:', '%.2e' % empS
    print 'lognormal prediction of S for birds:', '%.2e' % S2


    return



Fig3()
