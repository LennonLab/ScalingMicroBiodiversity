from __future__ import division
from __future__ import division
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


mydir = os.path.expanduser("~/GitHub/MicrobialScaling/")
mydir2 = os.path.expanduser("~/")



def Fig2(condition, ones, sampling):


    """ A figure demonstrating a strong abundance relationship across 30
    orders of magnitude in total abundance. The abundance of the most abundant
    species scales in a log-log fashion with the total abundance of the sample
    or system. """

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

    fs = 13 # font size used across figures
    Nlist, NmaxList, klist, datasets, radDATA = [[],[],[],[],[]]


    for name in os.listdir(mydir +'data/micro'):
        #if name in BadNames: continue
        if name in GoodNames: pass
        else: continue

        path = mydir+'data/micro/'+name+'/'+name+tail
        numlines = sum(1 for line in open(path))
        print name, numlines
        datasets.append([name, 'micro', numlines])

    for dataset in datasets:
        name, kind, numlines = dataset
        lines = []

        small = ['BIGN', 'BOVINE', 'CHU', 'LAUB', 'SED']
        big = ['HUMAN', 'CHINA', 'CATLIN', 'FUNGI', 'HYDRO']

        if name in small:
            lines = np.random.choice(range(1, numlines+1), 1000, replace=True)

        elif name in big:
            lines = np.random.choice(range(1, numlines+1), 2500, replace=True)

        elif name == 'TARA':
            lines = np.random.choice(range(1, numlines+1), 2500, replace=True)
        else:
            lines = np.random.choice(range(1, numlines+1), 2500, replace=True)


        path = mydir+'data/micro/'+name+'/'+name+tail

        for line in lines:
            data = linecache.getline(path, line)
            radDATA.append(data)

        klist.append('DarkCyan')

    for data in radDATA:

        data = data.split()
        name, kind, N, S, Var, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = data

        N = float(N)
        S = float(S)

        #if S < 10 or N < 11: continue # Min species richness

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

    #print f.summary()
    #print intercept, slope

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

    label1 = 'Dominance scaling law for microbial data compilation'
    label2 = 'Ranges of published $N_{max}$ and $N$'

    plt.fill_between(Nlist2, pred_ci_low, pred_ci_upp, color='r', lw=0.5, alpha=0.2)
    plt.text(2, 22, r'$N_{max}$'+ ' = '+str(round(10**intercept,2))+'*'+r'$N$'+'$^{'+str(round(slope,2))+'}$', fontsize=fs+4, color='Crimson', alpha=0.9)
    plt.text(2, 19,  r'$r^2$' + ' = ' +str("%.2f" % R2), fontsize=fs+4, color='0.2')
    plt.plot(X.tolist(), Y.tolist(), '--', c='red', lw=2, alpha=0.8, color='Crimson', label=label1)

    print 'r-squared and slope for RADs w/out inferred:', round(R2, 3), round(slope,3)


    plt.hexbin(Nlist, NmaxList, mincnt=1, gridsize = 50, bins='log', cmap=plt.cm.Reds_r)  #
    #plt.scatter(Nlist, NmaxList, color = 'LightCoral', alpha= 0.6 , s = 10, linewidths=0.5, edgecolor='Crimson')


    GO = np.log10([360.0*(10**26), 1010.0*(10**26)]) # estimated open ocean bacteria; Whitman et al. 1998
    Pm = np.log10([2.8*(10**27), 3.0*(10**27)]) # estimated Prochlorococcus; Flombaum et al. 2013
    Syn = np.log10([6.7*(10**26), 7.3*(10**26)]) # estimated Synechococcus; Flombaum et al. 2013

    Earth = np.log10([9.2*(10**29), 31.7*(10**29)]) # estimated bacteria on Earth; Kallmeyer et al. 2012
    SAR11 = np.log10([2.0*(10**28), 2.0*(10**28)]) # estimated percent abundance of SAR11; Morris et al. (2002)

    HGx = np.log10([0.5*(10**14), 1.5*(10**14)]) # estimated bacteria in Human gut; Berg (1996)
    HGy = np.log10([0.05*(10**min(HGx)), 0.15*(10**max(HGx))]) # estimated most abundant bacteria in Human gut; Turnbaugh et al. (2009), & Dethlefsen et al. (2008)

    COWx = np.log10([0.5*2.226*(10**15), 1.5*2.226*(10**15)]) # estimated bacteria in Cow rumen; LOW:   HIGH: Whitman et al. (1998)
    COWy = np.log10([0.09*(10**min(COWx)), .15*(10**max(COWx))]) # estimated dominance in Cow rumen; Stevenson and Weimer (2006)

    c = '0.2'
    ## EARTH
    x = [np.mean(Earth)]
    x_range = (max(Earth) - min(Earth))/2.0
    y = [np.mean([min(Pm), max(SAR11)])]
    y_range = (max(SAR11) - min(Pm))/2.0

    ax.text(8.5, max(SAR11)+0.2, r'$Prochlorococcus$ and Pelagibacterales', fontsize=fs+2, color = 'k')
    ax.text(max(Earth)+0.5, 26, 'Earth microbiome', fontsize=fs+2, color = 'k', rotation = 90)
    ax.axhline(y, 0, 0.90, ls = '--', c = '0.4')
    ax.axvline(x, 0, 0.85, ls = '--', c = '0.4')
    plt.errorbar(x, y, xerr=x_range, yerr=y_range, color='k', linewidth=1, label=label2)

    c = '0.4'
    ## GLOBAL OCEAN
    x = [np.mean(GO)]
    x_range = (max(GO) - min(GO))/2.0
    y = [np.mean(Pm)]
    y_range = (max(SAR11) - min(Pm))/2.0

    ax.text(7.5, min(Pm)-1.35, r'$Synechococcus$ and $Prochlorococcus$', fontsize=fs+2, color = 'k')
    ax.text(min(GO)-1, 22, 'Non-sediment ocean bacteria', fontsize=fs+2, color = 'k', rotation = 90)
    ax.axhline(y, 0, 0.85, ls = '--', c = '0.4')
    ax.axvline(x, 0, 0.83, ls = '--', c = '0.4')
    plt.errorbar(x, y, xerr=x_range, yerr=y_range, color='k', linewidth=1)

    ## HUMAN GUT
    x = [np.mean(HGx)]
    x_range = (max(HGx) - min(HGx))/2.0
    y = [np.mean(HGy)]
    y_range = (max(HGy) - min(HGy))/2.0

    ax.text(4, min(HGy)-1, 'Human gut', fontsize=fs+2, color = 'k')
    ax.text(min(HGx)-1, 8, 'Human gut', fontsize=fs+2, color = 'k', rotation = 90)
    ax.axhline(y, 0, 0.40, ls = '--', c = '0.4')
    ax.axvline(x, 0, 0.38, ls = '--', c = '0.4')
    plt.errorbar(x, y, xerr=x_range, yerr=y_range, color='k', linewidth=1)

    ## COW RUMEN
    x = [np.mean(COWx)]
    x_range = (max(COWx) - min(COWx))/2.0
    y = [np.mean(COWy)]
    y_range = (max(COWy) - min(COWy))/2.0

    ax.text(7, max(COWy)+0.3, '$Prevotella$', fontsize=fs+2, color = 'k')
    ax.text(max(COWx)+0.4, 11.2, 'Cow rumen', fontsize=fs+2, color = 'k', rotation = 90)
    ax.axhline(y, 0, 0.41, ls = '--', c = '0.4')
    ax.axvline(x, 0, 0.43, ls = '--', c = '0.4')
    plt.errorbar(x, y, xerr=x_range, yerr=y_range, color='k', linewidth=1)

    ax.text(5, -4.2, 'Number of reads or total abundance, '+ '$log$'+r'$_{10}$', fontsize=fs+4)
    ax.text(-2.5, 22, metric, fontsize=fs+4, rotation=90)

    plt.plot([0,32],[0,32], ls = '--', lw=2, c='0.7')
    #ax.text(18, 21, '1:1 line', fontsize=fs*1.0, rotation=40, color='0.7')

    plt.xlim(1, 33)
    plt.ylim(0, 32)

    plt.legend(bbox_to_anchor=(-0.015, 1, 1.025, .2), loc=10, ncol=1,
                                mode="expand",prop={'size':fs+1}, numpoints=1)

    if ones == False:
        plt.savefig(mydir+'/figs/Fig2/Locey_Lennon_2015_Fig2-'+condition+'_NoSingletons_'+str(sampling)+'.pdf', dpi=300, bbox_inches = "tight")
    if ones == True:
        plt.savefig(mydir+'/figs/Fig2/Locey_Lennon_2015_Fig2-'+condition+'_'+str(sampling)+'.pdf', dpi=300, bbox_inches = "tight")
        #plt.savefig(mydir+'/figs/Fig2/figure2-v2.pdf', dpi=300, bbox_inches = "tight")

    #plt.show()
    return



""" The following lines call figure functions to reproduce figures from the
    Locey and Lennon (2014) manuscript """

EMPcondition = ['closed']
Singletons = [True]
Samplings = [100]

for condition in EMPcondition:
    for ones in Singletons:
        for sampling in Samplings:
            Fig2(condition, ones, sampling)
