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

    BadNames = ['.DS_Store', 'EMPclosed', 'AGSOIL', 'SLUDGE', 'FECES', 'FUNGI']

    for name in os.listdir(mydir2 +'data/micro'):
        if name in BadNames: continue
        #if name == 'EMPopen': pass
        #else: continue

        path = mydir2+'data/micro/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt'
        #path = mydir2+'data/micro/'+name+'/'+name+'-SADMetricData.txt'

        numlines = sum(1 for line in open(path))
        print name, numlines
        datasets.append([name, 'micro', numlines])

    its = 1
    for i in range(its):
        print i
        for dataset in datasets:
            name, kind, numlines = dataset
            lines = []

            lines = []
            lines = np.random.choice(range(1, numlines+1), 2000, replace=True)

            #lines = random.sample(range(1, numlines+1), numlines)

            path = mydir2+'data/'+kind+'/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt'
            #path = mydir2+'data/'+kind+'/'+name+'/'+name+'-SADMetricData.txt'

            for line in lines:
                data = linecache.getline(path, line)
                radDATA.append(data)

            klist.append('DarkCyan')

    for data in radDATA:

        data = data.split()
        name, kind, N, S, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = data

        N = float(N)
        S = float(S)

        if S < 10 or N < 11: continue # Min species richness

        Nlist.append(float(np.log10(float(N))))
        NmaxList.append(float(np.log10(float(Nmax))))
        klist.append('DarkCyan')

    metric = 'Dominance, '+'log'+r'$_{10}$'

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

    print 'r-squared and slope for RADs w/out inferred:', round(R2, 3), round(slope,3)

    #label1 = 'EMP (heat mapped only): slope ='+str(round(slope,3))+', ' + r'$R^2$' + '=' +str(round(R2, 3))
    #ax.plot([0],[0], '-', c='Steelblue', lw=4, alpha=1, label=label1)

    plt.hexbin(Nlist, NmaxList, mincnt=1, gridsize = 80, bins='log', cmap=plt.cm.Reds_r)

    GO = log10(1110*10**26) # estimated open ocean bacteria; add reference
    Pm = log10(2.9*10**27) # estimated Prochlorococcus marinus; add reference
    Earth = log10(10**30) # estimated bacteria on Earth; add reference
    SAR11 = log10(2*10**28) # estimated Pelagibacter ubique; add reference
    Earth = log10(3.17 * 10**30) # estimated bacteria on Earth; add reference

    HGx = log10(10**14) # estimated bacteria in Human gut; add reference
    HGy = log10(0.1169*(10**14)) # estimated most abundant bacteria in Human gut; add reference
    # 0.0053
    COWx = log10(2.226*10**15) # estimated bacteria in Cow rumen; add reference
    COWy = log10((0.52/80)*(2.226*10**15)) # estimated dominance in Cow rumen; add reference
    #0.5/80

    #EMPx = log10(1252724704)
    #EMPy = log10(597974)

    c = '0.3'
    ax.text(11, SAR11+0.5, 'Ocean abundance of '+r'$Pelagibacter$'+' '+r'$ubique$', fontsize=fs+2, color = c)
    ax.axhline(SAR11, 0, 0.9, ls = '--', c = c)

    ax.text(11, Pm-1.15, 'Abundance of '+r'$Prochloroccus$', fontsize=fs+2, color = c)
    ax.axhline(Pm, 0, 0.88, ls = '--', c = c)

    ax.text(GO-1, 24, 'Non-sediment ocean bacteria', fontsize=fs+2, color = c, rotation = 90)
    ax.axvline(GO, 0, 0.86, ls = '--', c = c)

    ax.text(Earth+0.5, 26, 'Global abundance of bacteria', fontsize=fs+2, color = c, rotation = 90)
    ax.axvline(Earth, 0, 0.88, ls = '--', c = c)

    ax.text(2, HGy-1.5, 'Avg. in human guts', fontsize=fs+2, color = c)
    ax.axhline(HGy, 0, 0.40, ls = '--', c = c)

    ax.text(HGx-1, 8, 'Human gut', fontsize=fs+2, color = c, rotation = 90)
    ax.axvline(HGx, 0, 0.38, ls = '--', c = c)

    ax.text(4, COWy+0.6, 'Avg among '+r'$Prevotella$', fontsize=fs+2, color = c)
    ax.axhline(COWy, 0, 0.44, ls = '--', c = c)

    ax.text(COWx+0.4, 10.8, 'Cow rumen', fontsize=fs+2, color = c, rotation = 90)
    ax.axvline(COWx, 0, 0.41, ls = '--', c = c)

    ax.text(9, -3.17, 'Total abundance, '+ 'log'+r'$_{10}$', fontsize=fs*1.8)
    ax.text(-2.5, 22, metric, fontsize=fs*1.8, rotation=90)

    plt.scatter([GO], [Pm], color = '0.3', alpha= 1 , s = 40, linewidths=0.5, edgecolor='0.2')
    plt.scatter([Earth], [SAR11], color = '0.3', alpha= 1 , s = 40, linewidths=0.5, edgecolor='0.2')

    plt.scatter([HGx], [HGy], color = '0.3', alpha= 1 , s = 40, linewidths=0.5, edgecolor='0.2')
    plt.scatter([COWx], [COWy], color = '0.3', alpha= 1 , s = 40, linewidths=0.5, edgecolor='0.2')

    #plt.scatter([EMPx], [EMPy], color = 'b', alpha= 1 , s = 40, linewidths=0.5, edgecolor='0.2')
    #plt.scatter([COWx], [COWy], color = '0.3', alpha= 1 , s = 40, linewidths=0.5, edgecolor='0.2')


    plt.plot([0,32],[0,32], ls = '-', lw=2, c='0.7')
    ax.text(18, 21, '1:1 line', fontsize=fs*1.0, rotation=40, color='0.7')


    Nlist.extend([HGx, GO, Earth, COWx])
    NmaxList.extend([HGy, Pm, SAR11, COWy])

    """
    d = pd.DataFrame({'N': list(Nlist)})
    d['y'] = list(NmaxList)
    f = smf.ols('y ~ N', d).fit()

    R2 = f.rsquared
    pval = f.pvalues
    intercept = f.params[0]
    slope = f.params[1]

    print ': r-squared and slope for RADs with inferred:', round(R2,3), round(slope,3)
    """

    #label2 = 'EMP + inferred points: slope ='+str(round(slope,3))+', ' + r'$R^2$' + '=' +str(round(R2,3))

    #plt.legend(bbox_to_anchor=(-0.015, 1, 1.025, .2), loc=10, ncol=1,
    #                        mode="expand",prop={'size':fs+4})

    plt.xlim(1, 33)
    plt.ylim(0, 32)

    plt.savefig(mydir+'/figs/Locey_Lennon_2015_Fig2-OpenReference_NoSingletons.png', dpi=600, bbox_inches = "tight")
    plt.close()

    return




""" The following lines call figure functions to reproduce figures from the
    Locey and Lennon (2014) manuscript """

Fig2()
