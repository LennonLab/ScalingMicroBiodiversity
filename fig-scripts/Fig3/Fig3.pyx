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
from math import log10
import linecache



mydir = os.path.expanduser("~/GitHub/rare-bio/")
mydir2 = os.path.expanduser("~/GitHub/")





def BigS():

    """ A figure demonstrating a strong richness relationship across 10 or 11
    orders of magnitude in total abundance. Taxonomic richness of a sample
    scales in a log-log fashion with the total abundance of the sample.
    """

    fs = 10 # font size used across figures

    Nlist, Slist, klist, ESimplist, = [[], [], [], []]

    OrC = 'open' # is the microbial data (Earth Microbiome Project) going to
                    # represent closed or open reference OTU assignment
    ct = 0
    klist = []

    bigN = 0
    bigS = 0

    #radDATA = open(mydir+'output/EMP'+OrC+'-RADdata.txt','r')
    radDATA = open('/Users/lisalocey/Desktop/data/micro/EMP'+OrC+'/EMP'+OrC+'-SADMetricData_NoMicrobe1s.txt','r')
    for data in radDATA:

        data_list = data.split()
        #N, S, Evar, ESimp, ENee, EHeip, EQ, EPielou, BP, SimpDom, rareRel, rareOnes, skew = data_list
        name, kind, N, S, Evar, ESimp, ENee, EHeip, EQ, EPielou, BP, SimpDom, rareRel, rareOnes, skew = data_list

        Nlist.append(float(N))
        Slist.append(float(S))
        bigN += float(N)

        klist.append('DarkCyan')
        ct+=1

    metrics = [['Richness, '+'log'+r'$_{10}$', Slist]]

    bigN = np.log10(bigN)
    bigS = np.log10(5594412)

    print bigN, bigS

    fig = plt.figure()
    for index, i in enumerate(metrics):

        ax = fig.add_subplot(1, 1, 1)

        if index == 1: continue

        metric = i[0]
        metlist = i[1]

        RADListX = []
        RADListY = []

        rads = 0

        for j, k in enumerate(klist):

            if k == 'DarkCyan':
                rads += 1
                RADListX.append(Nlist[j])
                RADListY.append(metlist[j])

        # scatter plots
        indices = range(1000)
        random.shuffle(indices)

        RADListX = np.log10(RADListX)
        RADListY = np.log10(RADListY)

        RADListX = RADListX.tolist()
        RADListY = RADListY.tolist()

        Y = np.array(RADListY)
        X = np.array(RADListX)

        RADslope2, RADintercept2, RADrval2, RADpval2, RADstderr2 = sc.stats.linregress(X,Y)
        print metric,': r-squared and slope for RADs w/out inferred:',round(RADrval2**2,3), round(RADslope2,3)
        print RADslope2, RADintercept2, RADrval2, RADpval2, RADstderr2

        z = np.polyfit(X, Y, 1)
        p = np.poly1d(z)
        xp = np.linspace(1, 2, 1000)

        label1 = 'EMP sites (heat mapped): slope ='+str(round(RADslope2,3))+', ' + r'$R^2$' + '=' +str(round(RADrval2**2,3))
        ax.plot([0],[0], '-', c='Steelblue', lw=4, alpha=1, label=label1)

        plt.hexbin(RADListX, RADListY, mincnt=1, gridsize = 30, bins='log', cmap=plt.cm.Blues_r, label='EMP')

        # Adding in derived/inferred points
        ax.text(5, bigS * 1.02, 'Earth Microbiome Project', fontsize=fs+2, color = '0.4')
        ax.axhline(bigS, 0, bigN/10    , ls = '--', c = '0.4')

        #ax.text(bigN * 0.96, bigS-1.1, 'Samples of the EMP', fontsize=fs+2, color = '0.4', rotation = 90)
        ax.axvline(bigN, 0, bigS/8, ls = '--', c = '0.4')
        plt.scatter([bigN], [bigS], color = 'r', alpha= 1 , s = 40, linewidths=0.5, edgecolor='r')


        HMP_S = np.log10(27483)
        HMP_N = np.log10(22618041)

        ax.text(0.75, HMP_S * 1.02, 'Human Microbiome Project', fontsize=fs+2, color = '0.4')
        ax.axhline(HMP_S, 0, HMP_N/10    , ls = '--', c = '0.4')

        #ax.text(HMP_N * 0.96, HMP_S-1, 'Samples of the HMP', fontsize=fs+2, color = '0.4', rotation = 90)
        ax.axvline(HMP_N, 0, HMP_S/8, ls = '--', c = '0.4')
        plt.scatter([HMP_N], [HMP_S], color = 'r', alpha= 1 , s = 40, linewidths=0.5, edgecolor='r')

        ax.text(2, -1., 'Total abundance, '+ 'log'+r'$_{10}$', fontsize=fs*2)
        ax.text(-1, 6, 'OTU '+ metric, fontsize=fs*2, rotation=90)

        X = X.tolist()
        X.extend([bigN, HMP_N])
        Y = Y.tolist()
        Y.extend([bigS, HMP_S])

        RADslope, RADintercept, RADrval, RADpval, RADstderr = sc.stats.linregress(X, Y)
        print metric,': r-squared and slope for RADs w/ inferred:', round(RADrval**2,3), round(RADslope,3)

        z = np.polyfit(X, Y, 1)
        p = np.poly1d(z)
        xp = np.linspace(0, 32, 1000)

        label2 = 'EMP sites + reference pts: slope ='+str(round(RADslope,3))+', ' + r'$R^2$' + '=' +str(round(RADrval**2,3))
        ax.plot(xp, p(xp), '--', c='red', lw=2, label=label2)

        plt.legend(bbox_to_anchor=(-0.015, 1, 1.025, .2), loc=10, ncol=1,
                                mode="expand",prop={'size':fs+4})


        plt.xlim(0, 10)
        plt.ylim(0, 8)

    plt.savefig(mydir+'/figs/Locey_Lennon_2015_Fig4-'+OrC+'_REF.png', dpi=600,
                bbox_inches = "tight")
    plt.close()

    return




""" The following lines call figure functions to reproduce figures from the
    Locey and Lennon (2014) manuscript """

BigS()
