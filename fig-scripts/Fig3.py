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

mydir = os.path.expanduser("~/GitHub/MicrobialScaling/")
mydir2 = os.path.expanduser("~/")
pi = math.pi



def alpha2(a, N, Nmax, Nmin=1):
    y = sqrt(pi*Nmin*Nmax)/(2.0*a) * exp((a * log2(sqrt(Nmax/Nmin)))**2.0)
    y = y * exp((log(2.0)/(2.0*a))**2.0)
    y = y * erf(a * log2(sqrt(Nmax/Nmin)) - log(2.0)/(2.0*a))
    y += erf(a * log2(sqrt(Nmax/Nmin)) + log(2.0)/(2.0*a))
    y -= N

    return y # find alpha

def s2(a, Nmax, Nmin=1):
    return sqrt(pi)/a * exp( (a * log2(sqrt(Nmax/Nmin)))**2) # Using equation 10



def getNmax(N, b, slope):
    return 10 ** (b + slope*(log10(N)))

def expS(N, b, slope):
    return 10 ** (b + slope*(log10(N))) # 0.78 + 0.37*



def getS(Nrange, sb, sz, db, dz, guess, NmaxRange = [], predictNmax=True):

    Dlist = []
    Slist_ln = []
    Slist_SvN = []
    Nlist = []

    for i in range(1000):
        N = float(np.random.uniform(Nrange)[1])
        Nlist.append(N)

        Nmax = 0
        if predictNmax == True:
            Nmax = getNmax(N, db, dz)
        else:
            Nmax = np.random.uniform(NmaxRange)[1]

        Dlist.append(Nmax)
        Nmin = 1
        a = opt.fsolve(alpha2, guess, (N, Nmax, Nmin))[0]
        #print guess, a

        S2 = s2(a, Nmax, 1)
        Slist_ln.append(S2)

        S = expS(N, sb, sz)
        Slist_SvN.append(S)

    return [log10(Slist_ln), log10(Slist_SvN), log10(Dlist), log10(Nlist)]






def Fig3(condition, ones, sampling):

    """ A figure demonstrating a strong richness relationship across 10 or 11
    orders of magnitude in total abundance. Taxonomic richness of a sample
    scales in a log-log fashion with the total abundance of the sample.
    """

    fs = 10 # font size used across figures
    metric = 'Richness, '+'log'+r'$_{10}$'

    tail = str()
    if ones is False:
        tail = '-SADMetricData_NoMicrobe1s.txt'
    elif ones is True:
        tail = '-SADMetricData.txt'

    datasets = []
    GoodNames = []
    emp = str()

    if condition == 'open': emp = 'EMPopen'
    elif condition == 'closed': emp = 'EMPclosed'

    GoodNames = [emp, 'TARA', 'HMP', 'BIGN', 'BOVINE', 'CHU', 'LAUB', 'SED', 'HUMAN', 'CHINA', 'CATLIN', 'FUNGI']

    print '\n'

    its = 1
    d_blist = []
    d_zlist = []
    s_blist = []
    s_zlist = []

    for name in os.listdir(mydir +'data/micro'):
        if name in GoodNames: pass
        else: continue

        path = mydir+'data/micro/'+name+'/'+name+tail
        numlines = sum(1 for line in open(path))
        #print name, numlines
        datasets.append([name, 'micro', numlines])

    if sampling <= 500: its = 100
    else: its = 10

    for i in range(its):

        Nlist, Slist, klist, NmaxList = [[],[],[],[]]

        for dataset in datasets:
            radDATA = []
            name, kind, numlines = dataset
            lines = []

            small_mgrast = ['BIGN', 'BOVINE', 'CHU', 'LAUB', 'SED']
            big_mgrast = ['HUMAN', 'CHINA', 'CATLIN', 'FUNGI', 'HYDRO']

            if kind == 'micro':
                if name in small_mgrast:
                    lines = np.random.choice(range(1, numlines+1), 160, replace=True) # 40

                elif name in big_mgrast:
                    lines = np.random.choice(range(1, numlines+1), 400, replace=True) # 100

                else:
                    lines = np.random.choice(range(1, numlines+1), 400, replace=True) # 100

            path = mydir+'data/micro/'+name+'/'+name+tail

            for line in lines:
                data = linecache.getline(path, line)
                radDATA.append(data)

            ct = 0
            for data in radDATA:
                data = data.split()
                if data == []: continue
                name, kind, N, S, Var, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = data

                N = float(N)
                S = float(S)
                Nmax = float(Nmax)

                ct += 1
                Nlist.append(float(np.log10(N)))
                Slist.append(float(np.log10(S)))

                NmaxList.append(float(np.log10(Nmax)))
                klist.append('DarkCyan')

            #print name, ct


        Nlist, Slist, NmaxList = zip(*sorted(zip(Nlist, Slist, NmaxList)))
        Nlist = list(Nlist)
        Slist = list(Slist)
        NmaxList = list(NmaxList)

        # Regression for Dominance (Nmax) vs. N
        d = pd.DataFrame({'N': Nlist})
        d['Nmax'] = NmaxList
        f = smf.ols('Nmax ~ N', d).fit()

        R2 = f.rsquared
        pval = f.pvalues[0]
        intercept = f.params[0]
        slope = f.params[1]

        d_blist.append(intercept)
        d_zlist.append(slope)

        # Regression for Richness (S) vs. N
        d = pd.DataFrame({'N': Nlist})
        d['S'] = Slist
        f = smf.ols('S ~ N', d).fit()

        R2 = f.rsquared
        pval = f.pvalues[0]
        intercept = f.params[0]
        slope = f.params[1]

        s_blist.append(intercept)
        s_zlist.append(slope)

    sb = np.mean(s_blist)
    sz = np.mean(s_zlist)

    db = np.mean(d_blist)
    dz = np.mean(d_zlist)

    #print 'R2 for Nmax vs. N:', round(dR2, 3)
    print 'Nmax =', round(10**db, 2), '*', 'N^', round(dz, 2)
    #print 'R2 for S vs. N:', round(R2, 3)
    print 'S =', round(10**sb, 2), '*', 'N^', round(sz, 2),'\n'


    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    plt.text(2, 11.2, r'$S$'+ ' = '+str(round(10**sb, 1))+'*'+r'$N$'+'$^{'+str(round(sz, 2))+'},$' + r' $r^2$' + '=' +str(round(R2,2)), fontsize=fs+4, color='Crimson', alpha=0.9)

    # code for prediction intervals
    X = np.linspace(5, 32, 100)
    Y = f.predict(exog=dict(N=X))
    Nlist2 = Nlist + X.tolist()
    Slist2 = Slist + Y.tolist()

    d = pd.DataFrame({'N': list(Nlist2)})
    d['y'] = list(Slist2)
    f = smf.ols('y ~ N', d).fit()

    st, data, ss2 = summary_table(f, alpha=0.05)
    #fittedvalues = data[:,2]
    #pred_mean_se = data[:,3]
    pred_mean_ci_low, pred_mean_ci_upp = data[:,4:6].T
    pred_ci_low, pred_ci_upp = data[:,6:8].T

    plt.fill_between(Nlist2, pred_ci_low, pred_ci_upp, color='r', lw=0.5, alpha=0.2)
    z = np.polyfit(Nlist2, Slist2, 1)
    p = np.poly1d(z)
    xp = np.linspace(0, 32, 1000)

    label1 = 'Richness-abundance scaling relationship, $S$ = 7.6$N^{0.35}$'
    label2 = 'Predicted $S$ using the lognormal and published $N$ and $N_{max}$'
    label3 = 'Predicted $S$ using the lognormal, published $N$, and $N_{max}$ = 0.4$N^{0.93}$'

    plt.plot(xp, p(xp), '--', c='red', lw=2, alpha=0.8, color='Crimson', label=label1)
    plt.scatter(Nlist, Slist, color = 'LightCoral', alpha= 1 , s = 10, linewidths=0.5, edgecolor='Crimson')
    #plt.hexbin(Nlist, Slist, mincnt=1, gridsize = 80, bins='log', cmap=plt.cm.Reds_r, label='EMP')

    # Adding in derived/inferred points
    c = '0.3'

    GO = [3.6*(10**28), 10.1*(10**28)] # estimated open ocean bacteria; Whitman et al. 1998
    Pm = [2.8*(10**27), 3.0*(10**27)] # estimated Prochlorococcus; Flombaum et al. 2013
    Syn = [6.7*(10**26), 7.3*(10**26)] # estimated Synechococcus; Flombaum et al. 2013

    Earth = [9.2*(10**29), 31.7*(10**29)] # estimated bacteria on Earth; Kallmeyer et al. 2012
    SAR11 = [2.0*(10**28), 2.0*(10**28)] # estimated percent abundance of SAR11; Morris et al. (2002)

    HGx = [0.5*(10**14), 1.5*(10**14)] # estimated bacteria in Human gut; Berg (1996)
    HGy = [0.05*min(HGx), 0.15*max(HGx)] # estimated most abundant bacteria in Human gut; Turnbaugh et al. (2009), & Dethlefsen et al. (2008)

    COWx = [0.5*2.226*(10**15), 1.5*2.226*(10**15)] # estimated bacteria in Cow rumen; LOW:   HIGH: Whitman et al. (1998)
    COWy = [0.09*min(COWx), 0.15*max(COWx)] # estimated dominance in Cow rumen; Stevenson and Weimer (2006)


    ## PREDICTIONS OF S BASED ON THE EMPIRICAL S VS. N SCALING LAW, AND BASED
    ## ON THE LOGNORMAL PREDICTIVE FRAMEWORK OF CURTIS AND SLOAN USING
    ## 1.) THE ESTIMATED NMAX AND 2.) THE PREDICTED NMAX

    Ns = []
    Ss = []
    DomSs = []

    # Global Ocean estimates based on Whitman et al. (1998) and P. marinus (2012 paper)
    guess = 0.1019
    yrange = [min(Syn), max(Pm)]
    Slist_ln, Slist_SvN, Dlist, Nlist = getS(GO, sb, sz, db, dz, guess, yrange, predictNmax=False)

    S_ln = np.mean(Slist_ln)
    S_ln_sem = stats.sem(Slist_ln, ddof=1)
    S_SvN = np.mean(Slist_SvN)
    S_SvN_sem = stats.sem(Slist_SvN, ddof=1)
    Nmax = np.mean(Dlist)
    Nmax_sem = stats.sem(Dlist, ddof=1)
    avgN = np.mean(Nlist)
    avgN_sem = stats.sem(Nlist, ddof=1)
    Ss.append(S_ln)

    print 'scaling law prediction of S for Global Ocean:', '%.3e' % 10**(S_SvN)
    print 'lognormal prediction of S for Global Ocean, using estimated Nmax:', '%.3e' % 10**S_ln

    guess = 0.1019
    Slist_ln, Slist_SvN, Dlist, Nlist = getS(GO, sb, sz, db, dz, guess, yrange, predictNmax=True)

    S_ln = np.mean(Slist_ln)
    S_ln_sem = stats.sem(Slist_ln, ddof=1)
    S_SvN = np.mean(Slist_SvN)
    S_SvN_sem = stats.sem(Slist_SvN, ddof=1)
    Nmax = np.mean(Dlist)
    Nmax_sem = stats.sem(Dlist, ddof=1)
    avgN = np.mean(Nlist)
    avgN_sem = stats.sem(Nlist, ddof=1)

    print 'lognormal prediction of S for Global Ocean, using predicted Nmax:', '%.3e' % 10**S_ln
    #print 'P.m.:', '%.2e' % float(2.9*10**27), 'Nmax:', '%.2e' % 10**Nmax,'\n'

    S2 = float(S_ln)
    N = float(avgN)
    S_sem = float(4*S_ln_sem)
    N_sem = float(4*avgN_sem)

    ax.text(15, S2*0.93, 'Global Ocean', fontsize=fs+2, color = 'k')
    ax.axhline(S2, 0, 0.91, ls = '--', c = '0.6')
    ax.text(N-1, S2*.80, 'Global ocean', fontsize=fs+2, color = 'k', rotation = 90)
    ax.axvline(N, 0, 0.7, ls = '--', c = '0.6')
    #plt.scatter([N], [S2], color = '0.2', alpha= 1 , s = 60, linewidths=1, edgecolor='k')
    Ns.append(N)
    DomSs.append(S2)
    #plt.errorbar([N], [S2], xerr=N_sem, yerr=S_sem, color='k', linewidth=2)


    # Earth, i.e., Global estimates based on Kallmeyer et al. (2012) and SAR11 (2002 paper)
    guess = 0.1060
    yrange = [min(Pm), max(SAR11)]
    Slist_ln, Slist_SvN, Dlist, Nlist = getS(Earth, sb, sz, db, dz, guess, yrange, predictNmax=False)

    S_ln = np.mean(Slist_ln)
    S_ln_sem = stats.sem(Slist_ln, ddof=1)
    S_SvN = np.mean(Slist_SvN)
    S_SvN_sem = stats.sem(Slist_SvN, ddof=1)
    Nmax = np.mean(Dlist)
    Nmax_sem = stats.sem(Dlist, ddof=1)
    avgN = np.mean(Nlist)
    avgN_sem = stats.sem(Nlist, ddof=1)
    Ss.append(S_ln)

    #print 'average N and sem:' '%.3e' % 10**avgN, '%.3e' % 10**avgN_sem
    #print 'average Nmax and sem:' '%.3e' % 10**Nmax, '%.3e' % 10**Nmax_sem

    print '\nscaling law prediction of S for Earth:', '%.3e' % 10**S_SvN #,'%.3e' % 10**S_SvN_sem #, '%.3e' % S_SvN_CI
    print 'lognormal prediction of S for Earth, using estimated Nmax:', '%.3e' % 10**S_ln #, '%.3e' % 10**S_ln_sem#, '%.3e' % S_ln_CI

    guess = 0.1060
    Slist_ln, Slist_SvN, Dlist, Nlist = getS(Earth, sb, sz, db, dz, guess, yrange, predictNmax=True)

    S_ln = np.mean(Slist_ln)
    S_ln_sem = stats.sem(Slist_ln, ddof=1)
    S_SvN = np.mean(Slist_SvN)
    S_SvN_sem = stats.sem(Slist_SvN, ddof=1)
    Nmax = np.mean(Dlist)
    Nmax_sem = stats.sem(Dlist, ddof=1)
    avgN = np.mean(Nlist)
    avgN_sem = stats.sem(Nlist, ddof=1)

    print 'lognormal prediction of S for Earth, using predicted Nmax:', '%.3e' % 10**S_ln #, '%.3e' % 10**S_ln_sem#, '%.3e' % S_ln_CI
    #print 'SAR11:', '%.2e' % float(2.4*10**28), 'Nmax:', '%.2e' % Nmax,'\n'

    S2 = float(S_ln)
    N = float(avgN)
    S_sem = float(4*S_ln_sem)
    N_sem = float(4*avgN_sem)

    ax.text(25, S2*1.025, 'Earth', fontsize=fs+2, color = 'k')
    ax.axhline(S2, 0, 0.95, ls = '--', c = '0.6')
    ax.text(N-1, 8, 'Earth', fontsize=fs+2, color = 'k', rotation = 90)
    ax.axvline(N, 0, 0.85, ls = '--', c = '0.6')
    #plt.scatter([N], [S2], color = '0.2', alpha= 1 , s = 60, linewidths=1, edgecolor='k')
    #plt.errorbar([N], [S2], xerr=N_sem, yerr=S_sem, color='k', linewidth=2)
    Ns.append(N)
    DomSs.append(S2)


    # Human Gut
    guess = 0.1509
    Slist_ln, Slist_SvN, Dlist, Nlist = getS(HGx, sb, sz, db, dz, guess, HGy, predictNmax=False)

    S_ln = np.mean(Slist_ln)
    S_ln_sem = stats.sem(Slist_ln, ddof=1)
    S_SvN = np.mean(Slist_SvN)
    S_SvN_sem = stats.sem(Slist_SvN, ddof=1)
    Nmax = np.mean(Dlist)
    Nmax_sem = stats.sem(Dlist, ddof=1)
    avgN = np.mean(Nlist)
    avgN_sem = stats.sem(Nlist, ddof=1)
    Ss.append(S_ln)


    Slist_ln, Slist_SvN, Dlist, Nlist = getS(HGx, sb, sz, db, dz, guess, HGy, predictNmax=True)

    S_ln = np.mean(Slist_ln)
    S_ln_sem = stats.sem(Slist_ln, ddof=1)
    S_SvN = np.mean(Slist_SvN)
    S_SvN_sem = stats.sem(Slist_SvN, ddof=1)
    Nmax = np.mean(Dlist)
    Nmax_sem = stats.sem(Dlist, ddof=1)
    avgN = np.mean(Nlist)
    avgN_sem = stats.sem(Nlist, ddof=1)

    S2 = float(S_ln)
    N = float(avgN)
    S_sem = float(4*S_ln_sem)
    N_sem = float(4*avgN_sem)

    ax.text(3, S2*.9, 'Human Gut', fontsize=fs+2, color = 'k')
    ax.axhline(S2, 0, 0.41, ls = '--', c = '0.6')
    ax.text(N-1, 3.2, 'Human Gut', fontsize=fs+2, color = 'k', rotation = 90)
    ax.axvline(N, 0, 0.33, ls = '--', c = '0.6')
    #plt.scatter([N], [S2], color = '0.2', alpha= 1 , s = 60, linewidths=1, edgecolor='k')
    #plt.errorbar([N], [S2], xerr=N_sem, yerr=S_sem, color='k', linewidth=2)
    Ns.append(N)
    DomSs.append(S2)
    #print 'predS for Human Gut:', '%.3e' % 10**S2


    # Cow Rumen
    guess = 0.1
    Slist_ln, Slist_SvN, Dlist, Nlist = getS(COWx, sb, sz, db, dz, guess, COWy, predictNmax=False)

    S_ln = np.mean(Slist_ln)
    S_ln_sem = stats.sem(Slist_ln, ddof=1)
    S_SvN = np.mean(Slist_SvN)
    S_SvN_sem = stats.sem(Slist_SvN, ddof=1)
    Nmax = np.mean(Dlist)
    Nmax_sem = stats.sem(Dlist, ddof=1)
    avgN = np.mean(Nlist)
    avgN_sem = stats.sem(Nlist, ddof=1)
    Ss.append(S_ln)

    Slist_ln, Slist_SvN, Dlist, Nlist = getS(COWx, sb, sz, db, dz, guess, COWy, predictNmax=True)

    S_ln = np.mean(Slist_ln)
    S_ln_sem = stats.sem(Slist_ln, ddof=1)
    S_SvN = np.mean(Slist_SvN)
    S_SvN_sem = stats.sem(Slist_SvN, ddof=1)
    Nmax = np.mean(Dlist)
    Nmax_sem = stats.sem(Dlist, ddof=1)
    avgN = np.mean(Nlist)
    avgN_sem = stats.sem(Nlist, ddof=1)

    S2 = float(S_ln)
    N = float(avgN)
    S_sem = float(4*S_ln_sem)
    N_sem = float(4*avgN_sem)

    ax.text(6, S2*1.04, 'Cow Rumen', fontsize=fs+2, color = 'k')
    ax.axhline(S2, 0, 0.43, ls = '--', c = '0.6')
    ax.text(N+0.3, 4.2, 'Cow Rumen', fontsize=fs+2, color = 'k', rotation = 90)
    ax.axvline(N, 0, 0.38, ls = '--', c = '0.6')
    Ns.append(N)
    DomSs.append(S2)

    plt.scatter(Ns, Ss, color = '0.4', alpha= 1, s = 50, linewidths=2, edgecolor='k', label=label2)
    plt.scatter(Ns, DomSs, color = 'SkyBlue', alpha= 1, s = 50, linewidths=2, edgecolor='Steelblue', label=label3)
    #plt.errorbar([N], [S2], xerr=N_sem, yerr=S_sem, color='k', linewidth=2)

    ax.text(3, -1, 'Number of reads or total abundance, '+ '$log$'+r'$_{10}$', fontsize=fs*1.8)
    ax.text(-2.3, 12, 'OTU '+ metric, fontsize=fs*2, rotation=90)
    plt.xlim(1, 31)
    plt.ylim(0.8, 14)

    plt.legend(bbox_to_anchor=(-0.015, 1, 1.025, .2), loc=10, ncol=1,
                                mode="expand",prop={'size':fs+2})

    if ones == False:
        plt.savefig(mydir+'/figs/Fig3/Locey_Lennon_2015_Fig3-'+condition+'_NoSingletons_'+str(sampling)+'.png', dpi=600, bbox_inches = "tight")
    if ones == True:
        plt.savefig(mydir+'/figs/Fig3/Locey_Lennon_2015_Fig3-'+condition+'_'+str(sampling)+'.png', dpi=600, bbox_inches = "tight")

    #plt.show()

    return


""" The following lines call figure functions to reproduce figures from the
    Locey and Lennon (2014) manuscript """

EMPcondition = ['closed']
Singletons = [False]
Samplings = [1000]

for condition in EMPcondition:
    for ones in Singletons:
        for sampling in Samplings:
            Fig3(condition, ones, sampling)
