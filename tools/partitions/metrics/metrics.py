from __future__ import division
import sys
import numpy as np
from scipy.stats import gaussian_kde
from scipy import stats
import re
import math
import random, decimal


""" functions for calculating statistical metrics for vectors of integers """

def get_kdens_obs(partitions, metric):

    if metric == 'gini':
        D = get_kdens_obs_gini(partitions) # inequality
    
    elif metric == 'evar':
        D = get_kdens_obs_Evar(partitions) # Evar cannot be calculated when zeros = True
    
    elif metric == 'median':
        D = get_kdens_obs_MD(partitions)
    
    elif metric == 'variance':
        D = get_kdens_obs_var(partitions) # variance
    
    elif  metric == 'skewness':
        D = get_kdens_obs_skew(partitions) # skewness
    
    return D
    

def simplest_gini(x): #x is a vector of integers
    """This script was obtained from: https://subversion.american.edu/aisaac/notes/blinder_slides.xhtml.
    It yields Gini's coefficient of inequality, a common metric in economics for characterizing inequality
    in distributions of wealth"""
    #initialize yn and ysum
    yn, ysum, countx = 0.0, 0.0, 0
    #compute yN and ysum
    for xn in sorted(x):
      yn = (yn + xn)
      ysum = (ysum + yn)
      countx = (countx + 1)
    #compute area below Lorenz curve
    B = ysum / (countx * yn)
    return(1 - 2*B)
    
    
def gini_sample(SADs):
    """ Compute Gini's coefficient for each macrostate in a random sample """
    Gs = []
    for sad in SADs:
        G = simplest_gini(sad)
        Gs.append(G)
    return Gs
def get_kdens_obs_gini(partitions):
    """ Finds the kernel density function across a sample of partitions 
        for a given total (N) and number of parts (S) """
    ginis = gini_sample(partitions)
    density = gaussian_kde(ginis)
    n = len(ginis)
    xs = np.linspace(0.0,1.0,n)
    density.covariance_factor = lambda : .5
    density._compute_covariance()
    D = [xs,density(xs)]
    return D

def e_var(partition):
    """ Calculate Smith and Wilson's evenness index Evar """
    P = np.log(partition)
    S = len(partition)
    X = 0
    for x in P:
        X += (x - np.mean(P))**2/S
    evar = 1 - 2/math.pi*np.arctan(X) 
    return(evar)    


def vars_sample(partitions):
    """ Compute Evar for each partition in a random sample of partitions for N and S"""
    _vars = []
    for partition in partitions:
        _var = np.var(partition)
        _vars.append(_var)
    return _vars
def get_kdens_obs_var(partitions):
    """ Finds the kernel density function across a sample of partitions 
        for a given total (N) and number of parts (S) """
    _vars = vars_sample(partitions)
    density = gaussian_kde(_vars)
    n = len(_vars)
    _min = min(_vars)
    _max = max(_vars)
    xs = np.linspace(_min,_max,n)
    density.covariance_factor = lambda : .5
    density._compute_covariance()
    D = [xs,density(xs)]
    return D



def Evars_sample(partitions):
    """ Compute Evar for each partition in a random sample of partitions for N and S"""
    Evars = []
    for partition in partitions:
        Evar = e_var(partition)
        Evars.append(Evar)
    return Evars
def get_kdens_obs_Evar(partitions):
    """ Finds the kernel density function across a sample of partitions 
        for a given total (N) and number of parts (S) """
    Evars = Evars_sample(partitions)
    density = gaussian_kde(Evars)
    n = len(Evars)
    xs = np.linspace(0.0,1.0,n)
    density.covariance_factor = lambda : .5
    density._compute_covariance()
    D = [xs,density(xs)]
    return D

""" functions for statistical skewness """
def get_skews(partitions):
    """ Find the statistical skewnness for each integer partition in a sample """
    skews = []
    for partition in partitions:
        skews.append(stats.skew(partition))
    return skews
def get_kdens_obs_skew(partitions):
    """ Finds the kernel density function for the statistical skewnness across a sample of integer partitions for a given N and S """
    skews = get_skews(partitions)
    density = gaussian_kde(skews)
    n = len(skews)
    xs = np.linspace(float(min(skews)),float(max(skews)),n)
    density.covariance_factor = lambda : .5
    density._compute_covariance()
    D = [xs,density(xs)]
    return D

""" functions for examining median summand values across a sample of integer partitions """
def get_kdens_obs_MD(partitions):
    """ Finds the kernel density function for the median summmand across a sample of integer partitions for a given N and S """
    MDs = MDs_sample(partitions)
    density = gaussian_kde(MDs)
    n = len(MDs)
    xs = np.linspace(0.0,float(max(MDs)),n)
    density.covariance_factor = lambda : .5
    density._compute_covariance()
    D = [xs,density(xs)]
    return D
def MDs_sample(partitions):
    """ find the median summand value for each integer partition in a sample """
    MDs = []
    for partition in partitions:
        MDs.append(np.median(partition))
    return MDs
