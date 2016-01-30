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

mydir = os.path.expanduser("~/GitHub/MicrobialScaling/")


def Fig1(cutoffs, Ones):

    datasets = []
    GoodNames = ['LAUB', 'CHU', 'HYDRO', 'CATLIN']

    for cutoff in cutoffs:
        for name in GoodNames:
            if Ones == 'N':
                path = mydir+'data/micro/'+name+'/'+name+cutoff+'/'+name+cutoff+'-SADMetricData_NoMicrobe1s.txt'
            if Ones == 'Y':
	            path = mydir+'data/micro/'+name+'/'+name+cutoff+'/'+name+cutoff+'-SADMetricData.txt'

            num_lines = sum(1 for line in open(path))
            datasets.append([name, cutoff, 'micro', num_lines])
            print name, num_lines

    metrics = ['Rarity, '+r'$log_{10}$',
            'Dominance, '+r'$log_{10}$',
            'Evenness, ' +r'$log_{10}$',
            'Richness, ' +r'$log_{10}$']

    fig = plt.figure()
    for index, i in enumerate(metrics):

        metric = i
        fig.add_subplot(2, 2, index+1)
        fs = 10 # font size used across figures

        c97IntList, c97CoefList, c99IntList, c99CoefList, c95CoefList, c95IntList, R2List, metlist = [[], [], [], [], [], [], [], []]
        Nlist, Slist, Evarlist, ESimplist, klist, radDATA, BPlist, NmaxList, rareSkews, KindList, StdList = [[], [], [], [], [], [], [], [], [], [], []]
        #name, kind, N, S, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]

        its = 100
        for n in range(its):

            #name, kind, N, S, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
            Nlist, Slist, Evarlist, ESimplist, klist, radDATA, BPlist, NmaxList, rareSkews, KindList, StdList = [[], [], [], [], [], [], [], [], [], [], []]

            radDATA = []

            for dataset in datasets:

                name, cutoff, kind, numlines = dataset
                lines = []

                lines = np.random.choice(range(1, numlines+1), numlines, replace=False)

                if Ones == 'N':
                    path = mydir+'data/'+kind+'/'+name+'/'+name+cutoff+'/'+name+cutoff+'-SADMetricData_NoMicrobe1s.txt'
                if Ones == 'Y':
                    path = mydir+'data/'+kind+'/'+name+'/'+name+cutoff+'/'+name+cutoff+'-SADMetricData.txt'

                for line in lines:
                    data = linecache.getline(path, line)
                    dlist = cutoff+' '+data
                    radDATA.append(dlist)

            for data in radDATA:

                data = data.split()
                cutoff, name, kind, N, S, Var, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = data

                N = float(N)
                S = float(S)

                #if S < 10 or N < 11: continue

                Nlist.append(float(np.log10(N)))
                Slist.append(float(np.log10(S)))

                ESimplist.append(float(np.log10(float(ESimp))))
                KindList.append(cutoff)

                BPlist.append(float(BP))
                NmaxList.append(float(np.log10(float(Nmax))))

                # log-modulo transformation of skewnness
                lms = np.log10(np.abs(float(skew)) + 1)
                if skew < 0: lms = lms * -1
                rareSkews.append(float(lms))

                if cutoff == '99':
                    klist.append('Steelblue')

                elif cutoff == '97':
                    klist.append('Crimson')

                elif cutoff == '95':
                    klist.append('0.4')

            if index == 0: metlist = list(rareSkews)
            elif index == 1: metlist = list(NmaxList)
            elif index == 2: metlist = list(ESimplist)
            elif index == 3: metlist = list(Slist)

            # Multiple regression
            d = pd.DataFrame({'N': list(Nlist)})
            d['y'] = list(metlist)
            d['Kind'] = list(KindList)
            f = smf.ols('y ~ N * Kind', d).fit()

            print f.summary()
            #print f.params

            c95IntList.append(f.params[0])
            c95CoefList.append(f.params[3])

            if f.pvalues[1] < 0.05:  c97IntList.append(f.params[1] + f.params[0])
            else: c97IntList.append(f.params[0])

            if f.pvalues[4] < 0.05: c97CoefList.append(f.params[4] + f.params[3])
            else: c97CoefList.append(f.params[3])


            if f.pvalues[2] < 0.05:  c99IntList.append(f.params[2] + f.params[0])
            else: c99IntList.append(f.params[0])

            if f.pvalues[5] < 0.05: c99CoefList.append(f.params[5] + f.params[3])
            else: c99CoefList.append(f.params[3])

            R2List.append(f.rsquared)


        c95PIx, c95Fitted = [[],[]]
        c95CiH, c95CiL = [[],[]]

        c97PIx, c97Fitted = [[],[]]
        c97CiH, c97CiL = [[],[]]

        c99PIx, c99Fitted = [[],[]]
        c99CiH, c99CiL = [[],[]]

        c95ListX = []
        c95ListY = []
        c97ListX = []
        c97ListY = []
        c99ListX = []
        c99ListY = []

        for j, k in enumerate(KindList):
            if k == '99':
                c99ListX.append(Nlist[j])
                c99ListY.append(metlist[j])
            if k == '97':
                c97ListX.append(Nlist[j])
                c97ListY.append(metlist[j])
            if k == '95':
                c95ListX.append(Nlist[j])
                c95ListY.append(metlist[j])

        print metric
        lm = smf.ols('y ~ N * Kind', d).fit()
        print lm.summary()
        print '\n\n'

        st, data, ss2 = summary_table(lm, alpha=0.05)
        # ss2: Obs, Dep Var Population, Predicted Value, Std Error Mean Predict,
        # Mean ci 95% low, Mean ci 95% upp, Predict ci 95% low, Predict ci 95% upp,
        # Residual, Std Error Residual, Student Residual, Cook's D

        #fittedvalues = data[:,2]
        #predict_mean_se = data[:,3]
        predict_mean_ci_low, predict_mean_ci_upp = data[:,4:6].T
        predict_ci_low, predict_ci_upp = data[:,6:8].T


        for j, kval in enumerate(KindList):
            if kval == '99':
                c99CiH.append(predict_mean_ci_upp[j])
                c99CiL.append(predict_mean_ci_low[j])
                c99PIx.append(Nlist[j])
                c99Fitted.append(f.fittedvalues[j])
            if kval == '97':
                c97CiH.append(predict_mean_ci_upp[j])
                c97CiL.append(predict_mean_ci_low[j])
                c97PIx.append(Nlist[j])
                c97Fitted.append(f.fittedvalues[j])
            if kval == '95':
                c95CiH.append(predict_mean_ci_upp[j])
                c95CiL.append(predict_mean_ci_low[j])
                c95PIx.append(Nlist[j])
                c95Fitted.append(f.fittedvalues[j])


        c99PIx, c99Fitted, c99CiH, c99CiL, c97PIx, c97Fitted, c97CiH, c97CiL, c95PIx, c95Fitted, c95CiH, c95CiL = zip(*sorted(zip(c99PIx, c99Fitted, c99CiH, c99CiL, c97PIx, c97Fitted, c97CiH, c97CiL, c95PIx, c95Fitted, c95CiH, c95CiL)))

        plt.scatter(c99ListX, c99ListY, facecolor = 'none', alpha= 1 , s = 5, linewidths=0.5, edgecolor='Steelblue')
        plt.scatter(c97ListX, c97ListY, facecolor = 'none', alpha= 1 , s = 5, linewidths=0.5, edgecolor='Crimson')
        plt.scatter(c95ListX, c95ListY, facecolor = 'none', alpha= 1 , s = 5, linewidths=0.5, edgecolor='0.4')

        #plt.fill_between(c99PIx, c99CiL, c99CiH, color='b', lw=0.0, alpha=0.3)
        #plt.fill_between(c97PIx, c97CiL, c97CiH, color='r', lw=0.0, alpha=0.3)
        #plt.fill_between(c95PIx, c95CiL, c95CiH, color='k', lw=0.0, alpha=0.3)

        plt.plot(c99PIx, c99Fitted,  color='b', ls='--', lw=1, alpha=0.8)
        plt.plot(c97PIx, c97Fitted,  color='r', ls='--', lw=1, alpha=0.8)
        plt.plot(c95PIx, c95Fitted,  color='0.2', ls='--', lw=1, alpha=0.8)

        c99Int = round(np.mean(c99IntList), 2)
        c99Coef = round(np.mean(c99CoefList), 2)
        c97Int = round(np.mean(c97IntList), 2)
        c97Coef = round(np.mean(c97CoefList), 2)
        c95Int = round(np.mean(c95IntList), 2)
        c95Coef = round(np.mean(c95CoefList), 2)

        R2 = round(np.mean(R2List), 2)

        if index == 0:
            plt.ylim(-0.1, 2.0)
            plt.xlim(1, 6)
            plt.text(1.35, 1.7, r'$99%$'+ ' = '+str(round(10**c99Int,2))+'*'+r'$N$'+'$^{'+str(round(c99Coef,2))+'}$', fontsize=fs, color='Steelblue')
            plt.text(1.35, 1.5, r'$97%$'+ ' = '+str(round(10**c97Int,2))+'*'+r'$N$'+'$^{'+str(round(c97Coef,2))+'}$', fontsize=fs, color='Crimson')
            plt.text(1.35, 1.3, r'$95%$'+ ' = '+str(round(10**c95Int,2))+'*'+r'$N$'+'$^{'+str(round(c95Coef,2))+'}$', fontsize=fs, color='0.3')
            plt.text(1.35, 1.1,  r'$R^2$' + '=' +str(R2), fontsize=fs-1, color='k')

            plt.scatter([0],[-1], color = 'none', alpha = 1, s=10, linewidths=0.9, edgecolor='Steelblue', label= '99% (n='+str(len(c99ListY))+')')
            plt.scatter([0],[-1], color = 'none', alpha = 1, s=10, linewidths=0.9, edgecolor='Crimson', label= '97% (n='+str(len(c97ListY))+')')
            plt.scatter([0],[-1], color = 'none', alpha = 1, s=10, linewidths=0.9, edgecolor='0.3', label= '95% (n='+str(len(c95ListY))+')')
            plt.legend(bbox_to_anchor=(-0.04, 1.05, 2.48, .2), loc=10, ncol=3, mode="expand",prop={'size':fs+2})

        elif index == 1:

            plt.plot([0,7],[0,7], ls = '--', lw=1, c='0.7')
            #ax.text(18, 21, '1:1 line', fontsize=fs*1.0, rotation=40, color='0.7')
            plt.ylim(0, 6)
            plt.xlim(1, 6)

            plt.text(1.35, 5.1, r'$99%$'+ ' = '+str(round(10**c99Int,2))+'*'+r'$N$'+'$^{'+str(round(c99Coef,2))+'}$', fontsize=fs, color='Steelblue')
            plt.text(1.35, 4.6, r'$97%$'+ ' = '+str(round(10**c97Int,2))+'*'+r'$N$'+'$^{'+str(round(c97Coef,2))+'}$', fontsize=fs, color='Crimson')
            plt.text(1.35, 4.1, r'$95%$'+ ' = '+str(round(10**c95Int,2))+'*'+r'$N$'+'$^{'+str(round(c95Coef,2))+'}$', fontsize=fs, color='0.3')
            plt.text(1.35, 3.6,  r'$R^2$' + '=' +str(R2), fontsize=fs-1, color='k')

        elif index == 2:
            plt.ylim(-3.0, 0.0)
            plt.xlim(1, 6)

            plt.text(1.35, -1.8, r'$99%$'+ ' = '+str(round(10**c99Int,2))+'*'+r'$N$'+'$^{'+str(round(c99Coef,2))+'}$', fontsize=fs, color='Steelblue')
            plt.text(1.35, -2.1, r'$97%$'+ ' = '+str(round(10**c97Int,2))+'*'+r'$N$'+'$^{'+str(round(c97Coef,2))+'}$', fontsize=fs, color='Crimson')
            plt.text(1.35, -2.4, r'$95%$'+ ' = '+str(round(10**c95Int,2))+'*'+r'$N$'+'$^{'+str(round(c95Coef,2))+'}$', fontsize=fs, color='0.3')
            plt.text(1.35, -2.7,  r'$R^2$' + '=' +str(R2), fontsize=fs-1, color='k')


        elif index == 3:
            plt.ylim(0.9, 4.5)
            plt.xlim(1, 6)

            plt.text(1.35, 3.9, r'$99%$'+ ' = '+str(round(10**c99Int,2))+'*'+r'$N$'+'$^{'+str(round(c99Coef,2))+'}$', fontsize=fs, color='Steelblue')
            plt.text(1.35, 3.5, r'$97%$'+ ' = '+str(round(10**c97Int,2))+'*'+r'$N$'+'$^{'+str(round(c97Coef,2))+'}$', fontsize=fs, color='Crimson')
            plt.text(1.35, 3.1, r'$95%$'+ ' = '+str(round(10**c95Int,2))+'*'+r'$N$'+'$^{'+str(round(c95Coef,2))+'}$', fontsize=fs, color='0.3')
            plt.text(1.35, 2.7,  r'$R^2$' + '=' +str(R2), fontsize=fs-1, color='k')


        plt.xlabel('Number of reads, '+ '$log$'+r'$_{10}$', fontsize=fs)
        plt.ylabel(metric, fontsize=fs)
        plt.tick_params(axis='both', which='major', labelsize=fs-3)

    plt.subplots_adjust(wspace=0.4, hspace=0.4)

    if Ones =='N':
        plt.savefig(mydir+'figs/appendix/PercentCutoff/PercentCutoff_NoMicrobe1s.png', dpi=600, bbox_inches = "tight")
    elif Ones =='Y':
        plt.savefig(mydir+'figs/appendix/PercentCutoff/PercentCutoff.png', dpi=600, bbox_inches = "tight")

    #plt.show()
    plt.close()

    return


Fig1(cutoffs=['99','97','95'], Ones='Y')
#Fig1(cutoffs=['99','97','95'], Ones='N')
