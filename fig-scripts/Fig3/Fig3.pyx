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

import pandas as pd
#import patsy
import mpmath as mpm
from scipy.optimize import fsolve
from math import log10, erf
import linecache

mydir = os.path.expanduser("~/GitHub/rare-bio/")
mydir2 = os.path.expanduser("~/")


def Preston(N, Nmax):

    left = (2.0 * float(N))/(np.sqrt(np.pi) * float(Nmax))
    func = lambda a : left - (erf(np.log10(2)/a) / a)

    guess = 0.02 # alpha is often ~0.2
    a = fsolve(func, guess)

    expS = (np.sqrt(np.pi) / a) * np.exp( (np.log10(2.0)/(2.0*a))**2.0)
    return a[0], expS[0]



def get_EMP_SSADs():

    DATA = mydir2 + "data/micro/EMPclosed/EMPclosed-SSADdata.txt"

    SSADdict = {}

    with open(DATA) as f:

        for d in f:
            if d.strip():

                d = d.split()
                species = d[0]
                #sample = d[1]
                abundance = float(d[2])

                if abundance > 0:
                    if species in SSADdict:
                        SSADdict[species].append(abundance)
                    else:
                        SSADdict[species] = [abundance]



    SSADs = []
    SSADlist = SSADdict.items()

    S = len(SSADlist)
    N = 0
    for tup in SSADlist:

        SSAD = tup[1]
        if len(SSAD) >= 1:

            N += sum(SSAD)
            #SSAD.sort()
            #SSAD.reverse()
            #SSADs.append(SSAD)

    return [N, S]



def Fig3():

    """ A figure demonstrating a strong richness relationship across 10 or 11
    orders of magnitude in total abundance. Taxonomic richness of a sample
    scales in a log-log fashion with the total abundance of the sample.
    """

    fs = 10 # font size used across figures
    Nlist, Slist, klist, NmaxList, datasets, radDATA = [[],[],[],[],[],[]]
    metric = 'Richness, '+'log'+r'$_{10}$'

    BadNames = ['.DS_Store', 'EMPclosed', 'AGSOIL', 'SLUDGE', 'FECES', 'FUNGI']

    for name in os.listdir(mydir2 +'data/micro'):
        if name in BadNames: continue

        #if name == 'EMPopen' or name == 'HMP': pass
        #else: continue

        #path = mydir2+'data/micro/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt'
        path = mydir2+'data/micro/'+name+'/'+name+'-SADMetricData.txt'

        numlines = sum(1 for line in open(path))
        print name, numlines
        datasets.append([name, 'micro', numlines])

    its = 1
    for i in range(its):
        for dataset in datasets:
            name, kind, numlines = dataset

            #if name == 'EMPopen' or name == 'HMP': pass
            #else: continue

            lines = []
            lines = np.random.choice(range(1, numlines+1), 2000, replace=True)

            #path = mydir2+'data/'+kind+'/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt'
            path = mydir2+'data/'+kind+'/'+name+'/'+name+'-SADMetricData.txt'

            for line in lines:
                data = linecache.getline(path, line)
                radDATA.append(data)


    for data in radDATA:
        data = data.split()
        name, kind, N, S, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = data

        N = float(N)
        S = float(chao1)
        Nmax = float(Nmax)

        if S < 10 or N < 11: continue # Min species richness

        Nlist.append(float(np.log10(N)))
        Slist.append(float(np.log10(S)))
        NmaxList.append(float(np.log10(Nmax)))
        klist.append('DarkCyan')


    N_open_ones = 1315651204
    S_open_ones = 5594412
    N_open_noones = 1252725686
    S_open_noones = 2826534
    N_closed_ones = 654448644
    S_closed_ones = 69444
    N_closed_noones = 648525168
    S_closed_noones = 64658

    bigN = np.log10(N_open_noones)
    bigS = np.log10(S_open_noones)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    Nlist, Slist, NmaxList = zip(*sorted(zip(Nlist, Slist, NmaxList)))
    Nlist = list(Nlist)
    Slist = list(Slist)
    NmaxList = list(NmaxList)

    # Regression for Dominance (Nmax)
    d = pd.DataFrame({'N': Nlist})
    d['y'] = NmaxList
    f = smf.ols('y ~ N', d).fit()

    dR2 = f.rsquared
    dpval = f.pvalues
    dintercept = f.params[0]
    dslope = f.params[1]
    print 'intercept and slope for Nmax vs. N:', round(dintercept, 3), round(dslope,3)

    # Regression
    d = pd.DataFrame({'N': Nlist})
    d['y'] = Slist
    f = smf.ols('y ~ N', d).fit()

    R2 = f.rsquared
    pval = f.pvalues
    intercept = f.params[0]
    slope = f.params[1]

    print 'r-squared and slope for RADs w/out inferred:', round(R2, 3), round(slope,3)

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

    #label2 = 'Micro + reference pts: slope ='+str(round(slope,3))+', ' + r'$R^2$' + '=' +str(round(R2,3))
    plt.plot(xp, p(xp), '--', c='red', lw=2, alpha=0.8)#, label=label2)


    #label1 = 'sites (heat mapped): slope ='+str(round(slope, 3))+', ' + r'$R^2$' + '=' +str(round(R2, 3))
    #ax.plot([0],[0], '-', c='Steelblue', lw=4, alpha=1)#, label=label1)
    plt.hexbin(Nlist, Slist, mincnt=1, gridsize = 40, bins='log', cmap=plt.cm.Blues_r, label='EMP')


    # Adding in derived/inferred points

    # N and S for the Earth Microbiome Project
    plt.scatter([bigN], [bigS], color = 'orange', alpha= 1 , s = 40, linewidths=0.5, edgecolor='w', label='Earth Microbiome Project')

    # N and S for Earth,
    GO = float(1110*10**26) # estimated open ocean bacteria; add reference
    Pm = float(2.9*10**27) # estimated Prochlorococcus marinus; add reference
    Earth = float(9*10**29) # estimated bacteria on Earth; add reference
    SAR11 = float(2*10**28) # estimated Pelagibacter ubique; add reference

    # Global Ocean estimates based on Whitman et al. (1998) and dominance relationship
    N = float(GO)
    b =  -0.42
    z =  0.96
    Nmax = 10**(b + z * np.log10(N))
    alpha, S = Preston(N, Nmax)

    print 'predicted Nmax and S for the Global Ocean:', Nmax, S
    S = np.log10(S)
    N = np.log10(N)
    plt.scatter([N], [S], color = 'cyan', alpha= 1 , s = 40, linewidths=0.5, edgecolor='b', label='Global ocean; predicted Nmax = '+str(round(Nmax,2)))
    Nlist.append(N)
    Slist.append(S)


    # Global Ocean estimates based on Whitman et al. (1998) and P. marinus (2012 paper)
    N = float(GO)
    Nmax = float(Pm)
    alpha, S = Preston(N, Nmax)

    print 'predicted S for the Global Ocean:', S
    S = np.log10(S)
    N = np.log10(N)
    plt.scatter([N], [S], color = 'b', alpha= 1 , s = 40, linewidths=0.5, edgecolor='cyan', label='Global ocean; estimated Nmax = '+str(round(Nmax,2)))
    Nlist.append(N)
    Slist.append(S)

    # Global estimates based on Kallmeyer et al. (2012) and SAR11 (2002 paper)
    N = float(Earth)
    Nmax = float(SAR11)
    alpha, S = Preston(N, Nmax)

    print 'predicted S for Earth based on SAR11:', S
    #S = np.log10(S)
    #N = np.log10(N)
    plt.scatter([N], [S], color = 'DodgerBlue', alpha= 1 , s = 40, linewidths=1, edgecolor='green', label='Earth; estimated Nmax = '+str(Nmax))
    #Nlist.append(N)
    #Slist.append(S)


    # Global estimates based on Kallmeyer et al. (2012) and dominance relationship
    N = float(Earth)
    b = -0.42 # intercept
    z = 0.96 # slope
    Nmax = 10**(b + z * np.log10(N))
    alpha, S = Preston(N, Nmax)

    #print 'predicted Nmax for Earth:', Nmax

    print 'predicted S for Earth based on Nmax vs. N:', S
    S = np.log10(S)
    N = np.log10(N)
    plt.scatter([N], [S], color = 'green', alpha= 1 , s = 40, linewidths=1, edgecolor='DodgerBlue', label='Earth; predicted Nmax = '+str(Nmax))
    Nlist.append(N)
    Slist.append(S)

    # N and S for the Human Microbiome Project,
    S = np.log10(27483)
    N = np.log10(22618041)
    plt.scatter([N], [S], color = 'b', alpha= 1 , s = 40, linewidths=1, edgecolor='w', label='Human Microbiome Project')
    Nlist.append(N)
    Slist.append(S)

    # N and S for Hydrothermal vents (archaea), Brazelton et al. 2009
    S = np.log10(817)
    N = np.log10(167031)
    plt.scatter([N], [S], color = 'g', alpha= 1 , s = 40, linewidths=1, edgecolor='w', label='Sampled hydrothermal vents')
    Nlist.append(N)
    Slist.append(S)


    # N and S for sampled Cow Rumen, Jami & Mizrahi 2012
    S = np.log10(4896)
    N = np.log10(100000)
    plt.scatter([N], [S], color = 'm', alpha= 1 , s = 40, linewidths=1, edgecolor='w', label='Sampled cow rumen')
    Nlist.append(N)
    Slist.append(S)

    # Smapled N and S from Human Gut, from Eckburg et al. 2005
    N = np.log10(11831)
    S = np.log10(395)
    plt.scatter([N], [S], color = 'Lime', alpha= 1 , s = 40, linewidths=1, edgecolor='w', label='Sampled human gut')
    Nlist.append(N)
    Slist.append(S)

    ax.text(8, -2., 'Total abundance, '+ 'log'+r'$_{10}$', fontsize=fs*2)
    ax.text(-2.3, 18, 'OTU '+ metric, fontsize=fs*2, rotation=90)

    # Regression
    d = pd.DataFrame({'N': list(Nlist)})
    d['y'] = list(Slist)
    f = smf.ols('y ~ N', d).fit()

    R2 = f.rsquared
    pval = f.pvalues
    intercept = f.params[0]
    slope = f.params[1]

    print 'r-squared and slope for RADs with inferred:', round(R2, 3), round(slope,3)

    leg = plt.legend(loc=2, prop={'size':fs})
    leg.draw_frame(False)

    #plt.legend(bbox_to_anchor=(-0.015, 1, 1.025, .2), loc=10, ncol=1,
    #    mode="expand",prop={'size':fs+4}) # another option for making a legend

    plt.xlim(1, 31)
    plt.ylim(0.8, 25)

    plt.savefig(mydir+'/figs/Locey_Lennon_2015_Fig3-OpenReference_NoSingletons.png', dpi=600, bbox_inches = "tight")
    #plt.savefig(mydir+'/figs/Locey_Lennon_2015_Fig3-OpenReference.png', dpi=600, bbox_inches = "tight")
    #plt.savefig(mydir+'/figs/Locey_Lennon_2015_Fig3-ClosedReference_NoSingletons.png', dpi=600, bbox_inches = "tight")
    #plt.savefig(mydir+'/figs/Locey_Lennon_2015_Fig3-ClosedReference.png', dpi=600, bbox_inches = "tight")
    #plt.show()

    #plt.close()

    return




""" The following lines call figure functions to reproduce figures from the
    Locey and Lennon (2014) manuscript """

Fig3()
