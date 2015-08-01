from __future__ import division
#from __future__ import print_function
import  matplotlib.pyplot as plt

import numpy as np
import random
import scipy as sc
from scipy import stats
from scipy.optimize import curve_fit

import os
import sys
from scipy.stats.distributions import t

import statsmodels.stats.api as sms
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from statsmodels.stats.outliers_influence import summary_table

from math import log10
import itertools as it

import pandas as pd
import linecache


mydir = os.path.expanduser("~/GitHub/rare-bio/")
mydir2 = os.path.expanduser("~/")



def Fig2():

    """ A figure demonstrating a strong abundance relationship across 30
    orders of magnitude in total abundance. The abundance of the most abundant
    species scales in a log-log fashion with the total abundance of the sample
    or system. """

    fs = 10 # font size used across figures
    Nlist, NmaxList, klist, datasets, radDATA = [[],[],[],[],[]]

    #BadNames = ['.DS_Store', 'EMPclosed', 'AGSOIL', 'SLUDGE', 'FECES', 'FUNGI']
    GoodNames = ['MGRAST', 'HMP', 'EMPclosed']

    for name in os.listdir(mydir2 +'data/micro'):
        #if name in BadNames: continue
        if name in GoodNames: pass
        else: continue

        #path = mydir2+'data/micro/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt'
        path = mydir2+'data/micro/'+name+'/'+name+'-SADMetricData.txt'

        numlines = sum(1 for line in open(path))
        print name, numlines
        datasets.append([name, 'micro', numlines])

    for dataset in datasets:
        name, kind, numlines = dataset
        lines = []
        lines = np.random.choice(range(1, numlines+1), 10000, replace=True)
        #path = mydir2+'data/micro/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt'
        path = mydir2+'data/micro/'+name+'/'+name+'-SADMetricData.txt'

        for line in lines:
            data = linecache.getline(path, line)
            radDATA.append(data)

        klist.append('DarkCyan')

    for data in radDATA:

        data = data.split()
        name, kind, N, S, Var, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = data

        N = float(N)
        S = float(S)

        if S < 10 or N < 11: continue # Min species richness

        Nlist.append(float(np.log10(float(N))))
        NmaxList.append(float(np.log10(float(Nmax))))
        klist.append('DarkCyan')

    metric = 'Dominance, '+'$log$'+r'$_{10}$'

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    Nlist, NmaxList = zip(*sorted(zip(Nlist, NmaxList)))
    Nlist = list(Nlist)
    NmaxList = list(NmaxList)

    # Regression
    d = pd.DataFrame({'N': list(Nlist)})
    d['y'] = list(NmaxList)
    f = smf.ols('y ~ N', d).fit()

    R2 = f.rsquared
    pval = f.pvalues
    intercept = f.params[0]
    slope = f.params[1]

    print f.summary()
    print intercept, slope

    X = np.linspace(6, 40, 100)
    Y = f.predict(exog=dict(N=X))
    Nlist2 = Nlist + X.tolist()
    NmaxList2 = NmaxList + Y.tolist()

    d = pd.DataFrame({'N': list(Nlist2)})
    d['y'] = list(NmaxList2)
    f = smf.ols('y ~ N', d).fit()

    st, data, ss2 = summary_table(f, alpha=0.05)
    fittedvalues = data[:,2]
    pred_mean_se = data[:,3]
    pred_mean_ci_low, pred_mean_ci_upp = data[:,4:6].T
    pred_ci_low, pred_ci_upp = data[:,6:8].T

    plt.fill_between(Nlist2, pred_ci_low, pred_ci_upp, color='r', lw=0.5, alpha=0.2)
    plt.text(2, 22, r'$N_{max}$'+ ' = '+str(round(intercept,2))+'*'+r'$N$'+'$^{'+str(round(slope,2))+'}$', fontsize=fs+4, color='Crimson', alpha=0.9)
    plt.text(2, 20,  r'$R^2$' + '=' +str(round(R2,2)), fontsize=fs+4, color='0.2')

    print 'r-squared and slope for RADs w/out inferred:', round(R2, 3), round(slope,3)


    #label1 = 'EMP (heat mapped only): slope ='+str(round(slope,3))+', ' + r'$R^2$' + '=' +str(round(R2, 3))
    #ax.plot([0],[0], '-', c='Steelblue', lw=4, alpha=1, label=label1)

    plt.hexbin(Nlist, NmaxList, mincnt=1, gridsize = 80, bins='log', cmap=plt.cm.Reds_r)

    GO = np.log10([360.0*(10**26), 1010.0*(10**26)]) # estimated open ocean bacteria; Whitman et al. 1998
    Pm = np.log10([2.8*(10**27), 3.0*(10**27)]) # estimated Prochlorococcus; Flombaum et al. 2013
    Syn = np.log10([6.7*(10**26), 7.3*(10**26)]) # estimated Synechococcus; Flombaum et al. 2013

    Earth = np.log10([9.2*(10**29), 31.7*(10**29)]) # estimated bacteria on Earth; Kallmeyer et al. 2012
    SAR11 = np.log10([2.0*(10**28), 2.0*(10**28)]) # estimated percent abundance of SAR11; Morris et al. (2002)

    HGx = np.log10([0.5*(10**14), 1.5*(10**14)]) # estimated bacteria in Human gut; ADD REFERENCE
    HGy = np.log10([0.05*(10**min(HGx)), 0.15*(10**max(HGx))]) # estimated most abundant bacteria in Human gut; Turnbaugh et al. (2009)

    COWx = np.log10([0.5*2.226*(10**15), 1.5*2.226*(10**15)]) # estimated bacteria in Cow rumen; LOW:   HIGH: Whitman et al. (1998)
    COWy = np.log10([0.09*(10**min(COWx)), .15*(10**max(COWx))]) # estimated dominance in Cow rumen;

    c = '0.2'
    ## EARTH
    ax.text(11, max(SAR11+0.5), r'$Prochlorococcus$ to $Pelagibacterales$', fontsize=fs+2, color = c)
    ax.text(max(Earth)+0.5, 26, 'Earth microbiome', fontsize=fs+2, color = c, rotation = 90)
    ax.axhline(min(Pm), 0, 0.9, ls = '--', c = c)
    ax.axhline(max(SAR11), 0, 0.9, ls = '--', c = c)
    ax.axvline(min(Earth), 0, 0.88, ls = '--', c = c)
    ax.axvline(max(Earth), 0, 0.88, ls = '--', c = c)
    plt.fill_between(Earth, min(Pm), max(SAR11), color = c, alpha= 0.8 , linewidths=0.5, edgecolor='0.2')

    a = 0.8
    c = '0.4'
    ## GLOBAL OCEAN
    ax.text(10, min(Pm)-2.2, r'$Synechococcus$ to $Prochloroccus$', fontsize=fs+2, color = c)
    ax.text(min(GO)-1, 24, 'Non-sediment ocean bacteria', fontsize=fs+2, color = c, rotation = 90)
    ax.axhline(min(Syn), 0, 0.88, ls = '--', c = c)
    ax.axhline(max(Pm), 0, 0.88, ls = '--', c = c)
    ax.axvline(min(GO), 0, 0.86, ls = '--', c = c)
    ax.axvline(max(GO), 0, 0.86, ls = '--', c = c)
    plt.fill_between(GO, min(Syn), max(Pm), color = c, alpha= a , linewidths=0.5, edgecolor='0.2')

    c = '0.2'
    ## HUMAN GUT
    ax.text(2, min(HGy)-1.5, 'Human gut', fontsize=fs+2, color = c)
    ax.text(min(HGx)-1, 8, 'Human gut', fontsize=fs+2, color = c, rotation = 90)
    ax.axhline(min(HGy), 0, 0.40, ls = '--', c = c)
    ax.axhline(max(HGy), 0, 0.40, ls = '--', c = c)
    ax.axvline(min(HGx), 0, 0.38, ls = '--', c = c)
    ax.axvline(max(HGx), 0, 0.38, ls = '--', c = c)
    plt.fill_between(HGx, min(HGy), max(HGy), color = c, alpha= a, linewidths=0.5, edgecolor='0.2')

    c = '0.4'
    ## COW RUMEN
    ax.text(4, max(COWy)+0.6, 'Cow rumen', fontsize=fs+2, color = c)
    ax.text(max(COWx)+0.4, 10.8, 'Cow rumen', fontsize=fs+2, color = c, rotation = 90)
    ax.axhline(min(COWy), 0, 0.44, ls = '--', c = c)
    ax.axhline(max(COWy), 0, 0.44, ls = '--', c = c)
    ax.axvline(min(COWx), 0, 0.41, ls = '--', c = c)
    ax.axvline(max(COWx), 0, 0.41, ls = '--', c = c)
    plt.fill_between(COWx, min(COWy), max(COWy), color = c, alpha= a, linewidths=0.5, edgecolor='0.2')

    ax.text(3, -4.2, 'Number of reads or total abundance, '+ '$log$'+r'$_{10}$', fontsize=fs*1.8)
    ax.text(-2.5, 22, metric, fontsize=fs*1.8, rotation=90)

    #plt.scatter([EMPx], [EMPy], color = 'b', alpha= 1 , s = 40, linewidths=0.5, edgecolor='0.2')
    #plt.scatter([COWx], [COWy], color = '0.3', alpha= 1 , s = 40, linewidths=0.5, edgecolor='0.2')


    plt.plot([0,32],[0,32], ls = '-', lw=2, c='0.7')
    ax.text(18, 21, '1:1 line', fontsize=fs*1.0, rotation=40, color='0.7')

    plt.xlim(1, 33)
    plt.ylim(0, 32)

    #plt.savefig(mydir+'/figs/Fig2/Locey_Lennon_2015_Fig2-OpenReference_NoSingletons.png', dpi=600, bbox_inches = "tight")
    #plt.savefig(mydir+'/figs/Fig2/Locey_Lennon_2015_Fig2-ClosedReference_NoSingletons.png', dpi=600, bbox_inches = "tight")
    #plt.savefig(mydir+'/figs/Fig2/Locey_Lennon_2015_Fig2-OpenReference.png', dpi=600, bbox_inches = "tight")
    plt.savefig(mydir+'/figs/Fig2/Locey_Lennon_2015_Fig2-ClosedReference.png', dpi=600, bbox_inches = "tight")

    #plt.show()
    return



""" The following lines call figure functions to reproduce figures from the
    Locey and Lennon (2014) manuscript """

Fig2()
