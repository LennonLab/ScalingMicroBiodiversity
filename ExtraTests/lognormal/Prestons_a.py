from __future__ import division
#from bigfloat import BigFloat, sqrt, exp, log, log2, erf, const_pi
import numpy as np
import math
from numpy import log, log2, exp, sqrt,log10
from scipy.optimize import fsolve
import scipy.optimize as opt
import matplotlib.pyplot as plt
from scipy.special import erf
import sys

pi = math.pi

GO = [3.6*(10**28), 10.1*(10**28)] # estimated open ocean bacteria; Whitman et al. 1998
Pm = [2.8*(10**27), 3.0*(10**27)] # estimated Prochlorococcus; Flombaum et al. 2013
Syn = [6.7*(10**26), 7.3*(10**26)] # estimated Synechococcus; Flombaum et al. 2013

Earth = [9.2*(10**29), 31.7*(10**29)] # estimated bacteria on Earth; Kallmeyer et al. 2012
SAR11 = [2.0*(10**28), 2.0*(10**28)] # estimated percent abundance of SAR11; Morris et al. (2002)

HGx = 10**14 # estimated bacteria in Human gut; add reference
HGy = 0.1169*HGx # estimated most abundant bacteria in Human gut; add reference

AvianN = 2.82*10**11
AvianNmax = 3*10**9
AvianS = 10500



def alpha1(a, Nmax, Nt):
    return (sqrt(pi) * Nmax)/(2.0*a) * erf(log(2.0)/a) - Nt # find alpha

def s1(a):
    return sqrt(pi)/a * exp( (log(2.0)/(2.0*a))**2.0 ) # Using equation 8


def alpha2(a, N, Nmax, Nmin):
    y = sqrt(pi*Nmin*Nmax)/(2.0*a) * exp((a * log2(sqrt(Nmax/Nmin)))**2.0)
    y = y * exp((log(2.0)/(2.0*a))**2.0)
    y = y * erf(a * log2(sqrt(Nmax/Nmin)) - log(2.0)/(2.0*a)) + erf(a * log2(sqrt(Nmax/Nmin)) + log(2.0)/(2.0*a))
    y -= N

    return y # find alpha

def s2(a, Nmax, Nmin):
    return sqrt(pi)/a * exp( (a * log2(sqrt(Nmax/Nmin)))**2) # Using equation 10


def getNmax(N):
    return 10 ** (1.02*(log10(N)) - 0.71)

def empS(N, b=log10(3.92), slope=0.4): # macrobes: b = 0.86, slope = 0.23
    return 10 ** (b + slope*(log10(N)))


#N = float(AvianN)
#Nmax = AvianNmax
#Nmin = 1.0
#Nmax = getNmax(AvianN)

N = float(max(GO))
Nmax = float(max(Syn))
Nmin = 1.0
#Nmax = getNmax(N)

############################################### Assuming Nmin = 1

guess = 0.099
guess = 0.1019

a = opt.fsolve(alpha2, guess, (N, Nmax, Nmin))[0]
print guess, a

S2 = s2(a, Nmax, Nmin)
print 'S2:','%.3e' % S2 # predicted from lognormal

S = empS(N)
print 'empS:','%.3e' % S # predicted from scaling
