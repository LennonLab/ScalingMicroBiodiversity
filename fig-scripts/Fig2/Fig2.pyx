from __future__ import division
#from __future__ import print_function
import  matplotlib.pyplot as plt

import numpy as np
from scipy.stats import gaussian_kde
import random
import scipy as sc
from scipy import stats
from scipy.optimize import curve_fit

import os
import sys
from scipy.stats.distributions import t

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



def get_kdens(summands):
    """ Finds the kernel density function """

    density = gaussian_kde(summands)
    n = 1000 #len(summands)
    xs = np.linspace(float(min(summands)),float(max(summands)),n)
    density.covariance_factor = lambda : .4
    density._compute_covariance()
    D = [xs,density(xs)]
    return D




def Fig3():

    """ A figure demonstrating a strong abundance relationship across 30
    orders of magnitude in total abundance. The abundance of the most abundant
    species scales in a log-log fashion with the total abundance of the sample
    or system. """

    fs = 10 # font size used across figures
    #color = str()

    Nlist, Slist, Evarlist, ESimplist, ENeelist, EHeiplist, EQlist = [[], [], [],
                                                            [], [], [], []]
    klist, Shanlist, BPlist, SimpDomlist, SinglesList, tenlist, onelist = [[], [],
                                                        [], [], [], [], []]

    NmaxList, rareOnesList, rareRelList, rarePairList, rareSumOnesList = [[], [],
                                                                    [], [], []]

    OrC = 'open' # is the microbial data (Earth Microbiome Project) going to
                    # represent closed or open reference OTU assignment

    ct = 0
    klist = []

    radDATA = open(mydir+'output/EMP'+OrC+'-RADdata.txt','r')
    for data in radDATA:

        data_list = data.split()
        #N, S, ESimp, EHeip, BP, SimpDom, rareRel, rareOnes, rareSumOnes, Var, Evar = data_list
        N, S, Evar, ESimp, ENee, EHeip, EQ, EPielou, BP, SimpDom, rareRel, rareOnes, skew = data_list
        Nlist.append(float(N))
        Slist.append(float(S))

        Evarlist.append(float(Evar))
        NmaxList.append(float(BP)*float(N))

        klist.append('DarkCyan')

        ct+=1



    metrics = [['Dominance, '+'log'+r'$_{10}$', NmaxList]]


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
                if index == 0: RADListX.append(Nlist[j])
                else: RADListX.append(Slist[j])

                if metlist[j] == 0: RADListY.append(1)
                else: RADListY.append(metlist[j])


        # scatter plots
        indices = range(1000)
        random.shuffle(indices)

        RADListX = np.log10(RADListX)
        RADListY = np.log10(RADListY)

        RADListX = RADListX.tolist()
        RADListY = RADListY.tolist()


        #for i in indices:
        #    ax.scatter(RADListX[i], RADListY[i], color = 'SkyBlue', alpha= 0.5 , s = 40, linewidths=0.5, edgecolor='DarkCyan')

        Y = np.array(RADListY)
        X = np.array(RADListX)

        RADslope2, RADintercept2, RADrval2, RADpval2, RADstderr2 = sc.stats.linregress(X,Y)
        print metric,': r-squared and slope for RADs w/out inferred:',round(RADrval2**2,3), round(RADslope2,3)

        z = np.polyfit(X,Y,1)
        p = np.poly1d(z)
        xp = np.linspace(1, 2, 1000)

        label1 = 'EMP (heat mapped only): slope ='+str(round(RADslope2,3))+', ' + r'$R^2$' + '=' +str(round(RADrval2**2,3))
        ax.plot([0],[0], '-', c='Steelblue', lw=4, alpha=1, label=label1)

        plt.hexbin(RADListX, RADListY, mincnt=1, gridsize = 30, bins='log', cmap=plt.cm.Blues_r, label='EMP')

        # Adding in derived/inferred points
        GO = log10(1.2*10**29) # estimated open ocean bacteria
        Pm = log10(2.9*10**27) # estimated Prochlorococcus marinus

        Earth = log10(3.17 *10**30) # estimated bacteria on Earth
        SAR11 = log10(2 * 10**28) # estimated Pelagibacter ubique

        HGx = log10(10**13)
        HGy = log10(0.0053*(10**13))

        COWx = log10(2.226*10**15)
        COWy = log10((0.5/80)*(2.226*10**15))


        ax.text(11, SAR11+0.5, 'Ocean abundance of '+r'$Pelagibacter$'+' '+r'$ubique$', fontsize=fs+2, color = '0.4')
        ax.axhline(SAR11, 0, 0.88, ls = '--', c = '0.4')

        ax.text(11, Pm-1.1, 'Ocean abundance of '+r'$Prochloroccus$', fontsize=fs+2, color = '0.4')
        ax.axhline(Pm, 0, 0.86, ls = '--', c = '0.4')

        ax.text(GO-1, 24, 'Global abundance of ocean bacteria', fontsize=fs+2, color = '0.4', rotation = 90)
        ax.axvline(GO, 0, 0.86, ls = '--', c = '0.4')

        ax.text(Earth+0.5, 26, 'Global abundance of bacteria', fontsize=fs+2, color = '0.4', rotation = 90)
        ax.axvline(Earth, 0, 0.89, ls = '--', c = '0.4')

        ax.text(4.2, HGy+0.65, 'max for any OTU', fontsize=fs+2, color = '0.4')
        ax.axhline(HGy, 0, 0.38, ls = '--', c = '0.4')

        ax.text(HGx-1, 8, 'Human gut', fontsize=fs+2, color = '0.4', rotation = 90)
        ax.axvline(HGx, 0, 0.33, ls = '--', c = '0.4')

        ax.text(5, COWy+0.65, 'avg among '+r'$Prevotella$', fontsize=fs+2, color = '0.4')
        ax.axhline(COWy, 0, 0.42, ls = '--', c = '0.4')

        ax.text(COWx+0.4, 10.8, 'Cow rumen', fontsize=fs+2, color = '0.4', rotation = 90)
        ax.axvline(COWx, 0, 0.41, ls = '--', c = '0.4')

        ax.text(5.5, -3.15, 'Total abundance, '+ 'log'+r'$_{10}$', fontsize=fs*1.8)
        ax.text(-3.1, 24, metric, fontsize=fs*1.8, rotation=90)


        plt.scatter([GO], [Pm], color = 'r', alpha= 1 , s = 40, linewidths=0.5, edgecolor='r')
        plt.scatter([Earth], [SAR11], color = 'r', alpha= 1 , s = 40, linewidths=0.5, edgecolor='r')

        plt.scatter([HGx], [HGy], color = 'r', alpha= 1 , s = 40, linewidths=0.5, edgecolor='r')
        plt.scatter([COWx], [COWy], color = 'r', alpha= 1 , s = 40, linewidths=0.5, edgecolor='r')

        X = X.tolist()
        X.extend([HGx, GO, Earth, COWx])
        Y = Y.tolist()
        Y.extend([HGy, Pm, SAR11, COWy])

        RADslope, RADintercept, RADrval, RADpval, RADstderr = sc.stats.linregress(X,Y)
        print metric,': r-squared and slope for RADs w/ inferred:', round(RADrval**2,3), round(RADslope,3)

        z = np.polyfit(X,Y, 1)
        p = np.poly1d(z)
        xp = np.linspace(0, 32, 1000)

        label2 = 'EMP + inferred points: slope ='+str(round(RADslope,3))+', ' + r'$R^2$' + '=' +str(round(RADrval**2,3))
        ax.plot(xp, p(xp), '--', c='red', lw=2, label=label2)

        plt.legend(bbox_to_anchor=(-0.015, 1, 1.025, .2), loc=10, ncol=1,
                                mode="expand",prop={'size':fs+4})


        plt.xlim(0, 34)
        plt.ylim(0, 32)

    plt.savefig(mydir+'/figs/Locey_Lennon_2015_Fig3v2-'+OrC+'_REF.png', dpi=600,
                bbox_inches = "tight")
    plt.close()

    return




""" The following lines call figure functions to reproduce figures from the
    Locey and Lennon (2014) manuscript """

Fig3()
