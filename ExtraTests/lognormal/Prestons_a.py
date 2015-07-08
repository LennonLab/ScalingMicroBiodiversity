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

GO = 1110*10**26 # estimated open ocean bacteria; add reference
Pm = 2.9*10**27 # estimated Prochlorococcus marinus; add reference
#Earth = 3*10**30 # estimated bacteria on Earth; add reference
SAR11 = 2*10**28 # estimated Pelagibacter ubique; add reference
Earth = 3.17 * 10**30 # estimated bacteria on Earth; add reference

HGx = 10**14 # estimated bacteria in Human gut; add reference
HGy = 0.1169*HGx # estimated most abundant bacteria in Human gut; add reference

'''
Human Gut, V3, A	188071	3319	28984
Human Gut, V3, B	154446	2592	16400
Human Gut, V3, C	148182	5671	16785
Human Gut, V6, A	135570	1745	14440
Human Gut, V6, B	154698	1673	15342
Human Gut, V6, C	151626	1731	18530
EMP	1252724704	5594411	597974
HMP	22618041	27483	81814
'''


def alpha1(a):
    return (sqrt(pi) * Nmax)/(2.0*a) * erf(log(2.0)/a) - Nt # find alpha

def s1(a):
    return sqrt(pi)/a * exp( (log(2.0)/(2.0*a))**2.0 ) # Using equation 8


def alpha2(a):
    y = sqrt(pi*Nmin*Nmax)/(2.0*a) * exp((a * log2(sqrt(Nmax/Nmin)))**2.0)
    y = y * exp((log(2.0)/(2.0*a))**2.0)
    y = y * erf(a * log2(sqrt(Nmax/Nmin)) - log(2.0)/(2.0*a)) + erf(a * log2(sqrt(Nmax/Nmin)) + log(2.0)/(2.0*a))
    y -= Nt

    return y # find alpha

def s2(a):
    return sqrt(pi)/a * exp( (a * log2(sqrt(Nmax/Nmin)))**2) # Using equation 10


def getNmax(N):
    return 10 ** (0.5 + 0.93*(log10(N)))

def empS(N, b=0.639, slope=0.431):
    return 10 ** (0.639 + 0.431*(log10(N)))


Nt = float(Earth)
#Nmax = float(Pm)
Nmax = getNmax(Earth)
Nmin = 1.0

############################################### not assuming anything about Nmin
fig = plt.figure()
fig.add_subplot(2,2,1)

guess = 0.0005
a = opt.fsolve(alpha1, guess)[0]
#a = opt.newton(alpha1, guess)

S1 = s1(a)
#print '\nalpha1 =', a, 'f(alpha1) = ', '%.2e' % alpha1(a), 'S1:','%.3e' % S1

alist = np.linspace(0.99*a, a*1.2, 1000)
plt.plot(alist, alpha1(alist))
plt.xlabel("alpha1")
plt.ylabel("f(alpha1)")
plt.grid()

############################################### Assuming Nmin = 1
fig.add_subplot(2,2,2)

guess = 0.099
a = opt.fsolve(alpha2, guess)[0]
#a = opt.newton(alpha2, guess, maxiter=100)

S2 = s2(a)
print '\nalpha2 = ', a, 'f(alpha2) = ','%.2e' % alpha2(a), 'S2:','%.3e' % S2

S = empS(Nt)
print 'empS:','%.3e' % S


#alist = np.linspace(0.9999*a, a*1.0001, 1000)
#plt.plot(alist, alpha2(alist))
#plt.xlabel("alpha2")
#plt.ylabel("f(alpha2)")
#plt.grid()
#plt.show()
