from __future__ import division
import  matplotlib.pyplot as plt

import numpy as np
import random
import scipy as sc
from scipy import stats

import os
import sys
from scipy.stats.distributions import t

import statsmodels.stats.api as sms
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from statsmodels.stats.outliers_influence import summary_table

import itertools as it
import pandas as pd
from math import log10
import linecache

mydir = os.path.expanduser("~/GitHub/rare-bio/")
mydir2 = os.path.expanduser("~/")
sys.path.append(mydir2 + "tools/metrics")
import metrics

sys.path.append(mydir2 + "GitHub/partitions/partitions/Python/")
import partitions


def Fig1():

    metric_names = ['Rarity, '+r'$log_{10}$',
            'Dominance, '+r'$log_{10}$',
            'Evenness, ' +r'$log_{10}$',
            'Richness, ' +r'$log_{10}$']

    # Get some partitions...
    rads = []
    IN = mydir+'/output/partitions.txt'
    num_lines = sum(1 for line in open(IN))
    print 'lines:',num_lines
    for line in open(IN):
        rad = eval(line)
        #print 'rad:',rad
        rads.append(rad)

    fig = plt.figure()
    for index, i in enumerate(metric_names):

        metric = i
        print metric

        fig.add_subplot(2, 2, index+1)
        fs = 10 # font size used across figures

        Nlist, Slist, Evarlist, ESimplist, radDATA, BPlist, NmaxList, rareSkews, StdList = [[], [], [], [], [], [], [], [], []]

        for RAD in rads:

            N = sum(RAD)
            S = len(RAD)
            Nmax = max(RAD)

            if Nmax < 10 ** (0.9*(log10(N))): continue
            if S < 2 or N < 10: continue

            Nlist.append(np.log10(N))
            Slist.append(np.log10(S))

            # Evenness
            #Evar = metrics.e_var(RAD)
            ESimp = metrics.e_simpson(RAD)
            #EQ = metrics.EQ(RAD)
            #O = metrics.OE(RAD)
            #Camargo = 0.0 # metrics.camargo(RAD)   # Takes too long
            #ENee = metrics.NHC(RAD)
            #EPielou = metrics.e_pielou(RAD)
            #EHeip = metrics.e_heip(RAD)

            # Dominance
            BP = metrics.Berger_Parker(RAD)
            #SimpDom = metrics.simpsons_dom(RAD)

            #McN = metrics.McNaughton(RAD)

            # Rarity
            skew = stats.skew(RAD)
            #logskew = metrics.Rlogskew(RAD)

            # Preston's alpha and richness, from Curtis and Sloan (2002).
            # Estimating prokaryotic diversity and its limits. PNAS.
            #preston_a, preston_S = metrics.Preston(RAD)

            # Richness estimators
            #chao1, ace, jknife1, jknife2 = metrics.EstimateS1(RAD)
            #margalef = metrics.Margalef(RAD)
            #menhinick = metrics.Menhinick(RAD)

            ESimplist.append(np.log10(ESimp))
            BPlist.append(BP)
            NmaxList.append(np.log10(Nmax))

            # log-modulo transformation of skewnness
            lms = np.log10(np.abs(skew) + 1)
            if skew < 0: lms = lms * -1
            rareSkews.append(lms)

        if index == 0: metlist = list(rareSkews)
        elif index == 1: metlist = list(NmaxList)
        elif index == 2: metlist = list(ESimplist)
        elif index == 3: metlist = list(Slist)

        # OLS
        d = pd.DataFrame({'N': list(Nlist)})
        d['y'] = list(metlist)
        f = smf.ols('y ~ N', d).fit()
        print f.summary()

        Int = round(f.params[0], 2)
        Coef = round(f.params[1], 2)
        p = round(f.pvalues[0], 3)
        r2 = round(f.rsquared, 3)

        st, data, ss2 = summary_table(f, alpha=0.05)
        # ss2: Obs, Dep Var Population, Predicted Value, Std Error Mean Predict,
        # Mean ci 95% low, Mean ci 95% upp, Predict ci 95% low, Predict ci 95% upp,
        # Residual, Std Error Residual, Student Residual, Cook's D

        fittedvalues = data[:,2]
        predict_mean_se = data[:,3]
        predict_mean_ci_low, predict_mean_ci_upp = data[:,4:6].T
        predict_ci_low, predict_ci_upp = data[:,6:8].T

        Nlist, predict_mean_ci_low, predict_mean_ci_upp = zip(*sorted(zip(Nlist, predict_mean_ci_low, predict_mean_ci_upp)))

        plt.scatter(Nlist, metlist, color = 'LightCoral', alpha= 1 , s = 4, linewidths=0.5, edgecolor='Crimson')
        plt.fill_between(Nlist, predict_mean_ci_low, predict_mean_ci_upp, color='0.3', lw=0.0, alpha=0.6)

        if index == 0:
            plt.ylim(-0.1, 2.2)
            plt.xlim(1, 6)
            plt.text(1.35, 1.5, r'$E$'+ ' = '+str(Int)+'+'+str(Coef)+'*'+r'$N$', fontsize=fs-1, color='0.2')
            plt.text(1.35, 1.1,  r'$r^2$' + '=' +str(r2), fontsize=fs-1, color='k')

        elif index == 1:
            plt.ylim(0, 6)
            plt.xlim(1, 6)
            plt.text(1.35, 4.1, r'$D$'+ ' = '+str(round(Int,2))+'+'+str(round(Coef,2))+'*'+r'$N$', fontsize=fs-1, color='0.2')
            plt.text(1.35, 2.9,  r'$r^2$' + '=' +str(r2), fontsize=fs-1, color='k')

        elif index == 2:
            plt.ylim(-3., 0.0)
            plt.xlim(1, 6)
            plt.text(1.35, -2.4, r'$R$'+ ' = '+str(round(Int,2))+'+'+str(round(Coef,2))+'*'+r'$N$', fontsize=fs-1, color='0.2')
            plt.text(1.35, -2.7,  r'$r^2$' + '=' +str(r2), fontsize=fs-1, color='k')

        elif index == 3:
            plt.ylim(0.9, 4.5)
            plt.xlim(1, 6)
            plt.text(1.35, 3.1, r'$S$'+ ' = '+str(round(Int,2))+'+'+str(round(Coef,2))+'*'+r'$N$', fontsize=fs-1, color='0.2')
            plt.text(1.35, 2.3,  r'$r^2$' + '=' +str(r2), fontsize=fs-1, color='k')

        plt.xlabel('Total abundance, ' + r'$log_{10}$', fontsize=fs)
        plt.ylabel(metric, fontsize=fs)
        plt.tick_params(axis='both', which='major', labelsize=fs-3)

    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    plt.savefig(mydir+'figs/appendix/Constraints/Locey_Lennon_2015.png', dpi=600, bbox_inches = "tight")

    #plt.show()
    plt.close()
    return


Fig1()
