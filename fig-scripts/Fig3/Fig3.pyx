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
pi = math.pi


def alpha2(a, N, Nmax, Nmin):
    y = sqrt(pi*Nmin*Nmax)/(2.0*a) * exp((a * log2(sqrt(Nmax/Nmin)))**2.0)
    y = y * exp((log(2.0)/(2.0*a))**2.0)
    y = y * erf(a * log2(sqrt(Nmax/Nmin)) - log(2.0)/(2.0*a)) + erf(a * log2(sqrt(Nmax/Nmin)) + log(2.0)/(2.0*a))
    y -= N

    return y # find alpha

def s2(a, Nmax, Nmin):
    return sqrt(pi)/a * exp( (a * log2(sqrt(Nmax/Nmin)))**2) # Using equation 10



def getNmax(N):
    return 10 ** (-0.4 + 0.94*(log10(N)))


def expS(N, b, slope):
    return 10 ** (b + slope*(log10(N)))



def Fig3():

    """ A figure demonstrating a strong richness relationship across 10 or 11
    orders of magnitude in total abundance. Taxonomic richness of a sample
    scales in a log-log fashion with the total abundance of the sample.
    """

    fs = 10 # font size used across figures
    Nlist, Slist, klist, NmaxList, datasets, radDATA = [[],[],[],[],[],[]]
    metric = 'Richness, '+'log'+r'$_{10}$'

    #BadNames = ['.DS_Store', 'BCI', 'AGSOIL', 'SLUDGE', 'NABC', 'FECES', 'MGRAST', 'EMPopen']
    GoodNames = ['MGRAST', 'HMP', 'EMPopen']

    for name in os.listdir(mydir2 +'data/micro'):
        #if name in BadNames: continue
        if name in GoodNames: pass
        else: continue

        #path = mydir2+'data/micro/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt'
        path = mydir2+'data/micro/'+name+'/'+name+'-SADMetricData.txt'

        numlines = sum(1 for line in open(path))
        print name, numlines
        datasets.append([name, 'micro', numlines])

    print '\n'

    its = 1000
    for i in range(its):
        for dataset in datasets:

            name, kind, numlines = dataset
            lines = []
            if name == 'EMPclosed' or name == 'EMPopen':
                lines = np.random.choice(range(1, numlines+1), 100, replace=True)
            elif kind == 'micro': lines = np.random.choice(range(1, numlines+1), 100, replace=True)

            #path = mydir2+'data/'+kind+'/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt'
            path = mydir2+'data/'+kind+'/'+name+'/'+name+'-SADMetricData.txt'
            for line in lines:
                data = linecache.getline(path, line)
                radDATA.append(data)


    for data in radDATA:
        data = data.split()
        name, kind, N, S, Var, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = data

        N = float(N)
        S = float(S)
        #if S > 10**4: print name
        Nmax = float(Nmax)

        if S < 2 or N < 11: continue # Min species richness

        Nlist.append(float(np.log10(N)))
        Slist.append(float(np.log10(S)))

        NmaxList.append(float(np.log10(Nmax)))
        klist.append('DarkCyan')

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    Nlist, Slist, NmaxList = zip(*sorted(zip(Nlist, Slist, NmaxList)))
    Nlist = list(Nlist)
    Slist = list(Slist)
    NmaxList = list(NmaxList)

    # Regression for Dominance (Nmax) vs. N
    d = pd.DataFrame({'N': Nlist})
    d['Nmax'] = NmaxList
    f = smf.ols('Nmax ~ N', d).fit()

    dR2 = f.rsquared
    dpval = f.pvalues[0]
    dintercept = f.params[0]
    dslope = f.params[1]
    print 'R2 for Nmax vs. N:', round(dR2,3)

    # Regression for Richness (S) vs. N
    d = pd.DataFrame({'N': Nlist})
    d['S'] = Slist
    f = smf.ols('S ~ N', d).fit()

    R2 = f.rsquared
    pval = f.pvalues[0]
    intercept = f.params[0]
    slope = f.params[1]


    print 'R2 for S vs. N:', round(R2,3),'\n'
    plt.text(2, 10, r'$N_{max}$'+ ' = '+str(round(intercept,2))+'*'+r'$N$'+'$^{'+str(round(slope,2))+'}$', fontsize=fs+4, color='Crimson', alpha=0.9)
    plt.text(2, 9,  r'$R^2$' + '=' +str(round(R2,2)), fontsize=fs+4, color='0.2')

    # code for prediction intervals
    X = np.linspace(5, 32, 100)
    Y = f.predict(exog=dict(N=X))
    Nlist2 = Nlist + X.tolist()
    Slist2 = Slist + Y.tolist()

    d = pd.DataFrame({'N': list(Nlist2)})
    d['y'] = list(Slist2)
    f = smf.ols('y ~ N', d).fit()

    st, data, ss2 = summary_table(f, alpha=0.05)
    fittedvalues = data[:,2]
    pred_mean_se = data[:,3]
    pred_mean_ci_low, pred_mean_ci_upp = data[:,4:6].T
    pred_ci_low, pred_ci_upp = data[:,6:8].T

    plt.fill_between(Nlist2, pred_ci_low, pred_ci_upp, color='r', lw=0.5, alpha=0.2)
    z = np.polyfit(Nlist2, Slist2, 1)
    p = np.poly1d(z)
    xp = np.linspace(0, 32, 1000)

    plt.plot(xp, p(xp), '--', c='red', lw=2, alpha=0.8, label= r'$S$'+ ' = '+str(round(intercept,2))+'+'+str(round(slope,2))+'*'+r'$N$', color='Crimson')
    plt.hexbin(Nlist, Slist, mincnt=1, gridsize = 20, bins='log', cmap=plt.cm.Reds, label='EMP')

    # Adding in derived/inferred points
    c = '0.3'
    GO = 1110*10**26 # estimated open ocean bacteria; add reference
    Pm = 2.9*10**27 # estimated Prochlorococcus marinus; add reference
    Earth = 3.17*10**30 # estimated bacteria on Earth; add reference
    SAR11 = 2*10**28 # estimated Pelagibacter ubique; add reference

    HGx =10**14 # estimated bacteria in Human gut; add reference
    HGy = 0.1169*(10**14) # estimated most abundant bacteria in Human gut; add reference # 0.0053
    COWx = 2.226*10**15 # estimated bacteria in Cow rumen; add reference
    COWy = (0.52/80)*(2.226*10**15) # estimated dominance in Cow rumen; add reference #0.5/80

    # Global Ocean estimates based on Whitman et al. (1998) and P. marinus (2012 paper)
    Nmin = 1
    N = float(GO)
    empS = expS(N, intercept, slope)

    Nmax = Pm
    guess = 0.1019
    a = opt.fsolve(alpha2, guess, (N, Nmax, Nmin))[0]
    print guess, a
    S2 = s2(a, Nmax, Nmin)

    print 'P.m.:', '%.2e' % float(Pm), 'Nmax:', '%.2e' % getNmax(N)
    print 'scaling law prediction of S for Global Ocean:', '%.3e' % empS
    print 'lognormal prediction of S for Global Ocean, using estimated Nmax:', '%.3e' % S2

    Nmax = getNmax(N)
    guess = 0.106
    a = opt.fsolve(alpha2, guess, (N, Nmax, Nmin))[0]
    S2 = s2(a, Nmax, Nmin)

    print 'lognormal prediction of S for Global Ocean, using predicted Nmax:', '%.3e' % S2,'\n'

    S2 = log10(S2)
    N = log10(N)

    ax.text(17, S2*0.95, 'Global Ocean', fontsize=fs+2, color = c)
    ax.axhline(S2, 0, 0.93, ls = '--', c = c)
    ax.text(N-1, S2*.75, 'Global ocean', fontsize=fs+2, color = c, rotation = 90)
    ax.axvline(N, 0, 0.8, ls = '--', c = c)
    plt.scatter([N], [S2], color = '0.2', alpha= 1 , s = 60, linewidths=1, edgecolor='k')
    Nlist.extend([N])
    Slist.extend([S2])


    # Global estimates based on Kallmeyer et al. (2012) and SAR11 (2002 paper)
    N = float(Earth)
    empS = expS(N, intercept, slope)
    print slope, intercept, '%.3e' % empS

    Nmax = SAR11
    guess = 0.1060
    a = opt.fsolve(alpha2, guess, (N, Nmax, Nmin))[0]
    print guess, a
    S2 = s2(a, Nmax, Nmin)

    print 'P.ubique.:', '%.2e' % float(SAR11), 'Nmax:', '%.2e' % getNmax(N)
    print 'scaling law prediction of S for Earth:', '%.3e' % empS
    print 'lognormal prediction of S for Earth, using estimated Nmax:', '%.3e' % S2


    Nmax = getNmax(N)
    guess = 0.1011
    a = opt.fsolve(alpha2, guess, (N, Nmax, Nmin))[0]
    S2 = s2(a, Nmax, Nmin)

    print 'lognormal prediction of S for Earth, using predicted Nmax:', '%.3e' % S2

    S2 = log10(S2)
    N = log10(N)

    ax.text(20, S2*1.025, 'Earth', fontsize=fs+2, color = c)
    ax.axhline(S2, 0, 0.97, ls = '--', c = c)
    ax.text(N-1, 8, 'Earth', fontsize=fs+2, color = c, rotation = 90)
    ax.axvline(N, 0, 0.8, ls = '--', c = c)
    plt.scatter([N], [S2], color = '0.2', alpha= 1 , s = 60, linewidths=1, edgecolor='k')
    Nlist.extend([N])
    Slist.extend([S2])


    # Human Gut based on ...
    N = float(HGx)
    #Nmax = float(HGy)
    Nmax = getNmax(N)

    guess = 0.1509
    a = opt.fsolve(alpha2, guess, (N, Nmax, Nmin))[0]
    S2 = s2(a, Nmax, Nmin)
    S2 = log10(S2)
    N = log10(N)

    ax.text(2, S2*.9, 'Human Gut', fontsize=fs+2, color = c)
    ax.axhline(S2, 0, 0.43, ls = '--', c = c)
    ax.text(N-1, 3.2, 'Human Gut', fontsize=fs+2, color = c, rotation = 90)
    ax.axvline(N, 0, 0.33, ls = '--', c = c)
    plt.scatter([N], [S2], color = '0.2', alpha= 1 , s = 60, linewidths=1, edgecolor='k')
    Nlist.extend([N])
    Slist.extend([S2])
    #print 'predS for Human Gut:', '%.3e' % 10**S2


    # Cow Rumen based on ...
    N = float(COWx)
    #Nmax = float(COWy)
    Nmax = getNmax(N)

    guess = 0.1
    a = opt.fsolve(alpha2, guess, (N, Nmax, Nmin))[0]
    S2 = s2(a, Nmax, Nmin)
    S2 = log10(S2)
    N = log10(N)

    ax.text(2, S2*1.04, 'Cow Rumen', fontsize=fs+2, color = c)
    ax.axhline(S2, 0, 0.45, ls = '--', c = c)
    ax.text(N+0.3, 4.2, 'Cow Rumen', fontsize=fs+2, color = c, rotation = 90)
    ax.axvline(N, 0, 0.38, ls = '--', c = c)
    plt.scatter([N], [S2], color = '0.2', alpha= 1 , s = 60, linewidths=1, edgecolor='k')
    Nlist.extend([N])
    Slist.extend([S2])


    ax.text(3, -1, 'Number of reads or total abundance, '+ '$log$'+r'$_{10}$', fontsize=fs*1.8)
    ax.text(-2.3, 12, 'OTU '+ metric, fontsize=fs*2, rotation=90)
    #leg = plt.legend(loc=2, numpoints = 1, prop={'size':fs})
    #leg.draw_frame(False)
    plt.xlim(1, 31)
    plt.ylim(0.8, 14)

    #plt.savefig(mydir+'/figs/Fig3/Locey_Lennon_2015_Fig3_OpenReference_NoSingletons.png', dpi=600, bbox_inches = "tight")
    plt.savefig(mydir+'/figs/Fig3/Locey_Lennon_2015_Fig3-OpenReference.png', dpi=600, bbox_inches = "tight")
    #plt.savefig(mydir+'/figs/Fig3/Locey_Lennon_2015_Fig3-ClosedReference_NoSingletons.png', dpi=600, bbox_inches = "tight")
    #plt.savefig(mydir+'/figs/Fig3/Locey_Lennon_2015_Fig3-ClosedReference.png', dpi=600, bbox_inches = "tight")
    #plt.show()

    return


""" The following lines call figure functions to reproduce figures from the
    Locey and Lennon (2015) manuscript """

Fig3()
