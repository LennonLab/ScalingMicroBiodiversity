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



def Fig1():

    datasets = []
    BadNames = ['.DS_Store', 'BCI', 'AGSOIL', 'SLUDGE', 'NABC']

    for name in os.listdir(mydir2 +'data/micro'):
        if name in BadNames: continue

        #path = mydir2+'data/micro/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt'
        path = mydir2+'data/micro/'+name+'/'+name+'-SADMetricData.txt'

        num_lines = sum(1 for line in open(path))
        datasets.append([name, 'micro', num_lines])
        print name, num_lines

    for name in os.listdir(mydir2 +'data/macro'):
        if name in BadNames: continue

        #path = mydir2+'data/macro/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt'
        path = mydir2+'data/macro/'+name+'/'+name+'-SADMetricData.txt'

        num_lines = sum(1 for line in open(path))
        datasets.append([name, 'macro', num_lines])
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

        MicIntList, MicCoefList, MacIntList, MacCoefList, R2List, metlist = [[], [], [], [], [], []]
        Nlist, Slist, Evarlist, ESimplist, klist, radDATA, BPlist, NmaxList, rareSkews, KindList, StdList = [[], [], [], [], [], [], [], [], [], [], []]
        #name, kind, N, S, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]

        its = 300
        for n in range(its):

            #name, kind, N, S, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
            Nlist, Slist, Evarlist, ESimplist, klist, radDATA, BPlist, NmaxList, rareSkews, KindList, StdList = [[], [], [], [], [], [], [], [], [], [], []]

            numMac = 0
            numMic = 0
            radDATA = []

            for dataset in datasets:

                name, kind, numlines = dataset
                lines = []
                if name == 'EMPclosed' or name == 'EMPopen':
                    lines = np.random.choice(range(1, numlines+1), 25, replace=True)
                if kind == 'micro': lines = np.random.choice(range(1, numlines+1), 50, replace=True)
                else: lines = np.random.choice(range(1, numlines+1), 80, replace=True)

                #path = mydir2+'data/'+kind+'/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt'
                path = mydir2+'data/'+kind+'/'+name+'/'+name+'-SADMetricData.txt'

                for line in lines:
                    data = linecache.getline(path, line)
                    radDATA.append(data)

            for data in radDATA:

                data = data.split()
                name, kind, N, S, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = data

                N = float(N)
                S = float(S)

                if S < 2 or N < 10: continue

                Nlist.append(float(np.log10(N)))
                Slist.append(float(np.log10(S)))

                ESimplist.append(float(np.log10(float(ESimp))))
                KindList.append(kind)

                BPlist.append(float(BP))
                NmaxList.append(float(np.log10(float(Nmax))))

                # log-modulo transformation of skewnness
                lms = np.log10(np.abs(float(skew)) + 1)
                if skew < 0: lms = lms * -1
                rareSkews.append(float(lms))

                if kind == 'micro':
                    numMic += 1
                    klist.append('b')
                if kind == 'macro':
                    klist.append('r')
                    numMac += 1

            if index == 0: metlist = list(rareSkews)
            elif index == 1: metlist = list(NmaxList)
            elif index == 2: metlist = list(ESimplist)
            elif index == 3: metlist = list(Slist)

            # Multiple regression
            d = pd.DataFrame({'N': list(Nlist)})
            d['y'] = list(metlist)
            d['Kind'] = list(KindList)
            f = smf.ols('y ~ N * Kind', d).fit()

            # Simple regression
            #f1 = smf.ols('y ~ N', d).fit()

            MacIntList.append(f.params[0])
            MacCoefList.append(f.params[2])

            if f.pvalues[1] < 0.05:
                MicIntList.append(f.params[1] + f.params[0])
            else:
                MicIntList.append(f.params[0])

            if f.pvalues[3] < 0.05:
                MicCoefList.append(f.params[3] + f.params[2])
            else:
                MicCoefList.append(f.params[2])

            R2List.append(f.rsquared)


        MacPIx, MacFitted, MicPIx, MicFitted = [[],[],[],[]]
        macCiH, macCiL, micCiH, micCiL = [[],[],[],[]]

        MacListX = []
        MacListY = []
        MicListX = []
        MicListY = []

        for j, k in enumerate(KindList):
            if k == 'micro':
                MicListX.append(Nlist[j])
                MicListY.append(metlist[j])

            elif k == 'macro':
                MacListX.append(Nlist[j])
                MacListY.append(metlist[j])

        print metric
        lm = smf.ols('y ~ N * Kind', d).fit()
        print lm.summary()
        print metric
        f1 = smf.ols('y ~ N', d).fit()
        print f1.summary()
        print '\n\n'

        st, data, ss2 = summary_table(lm, alpha=0.05)
        st_all, data_all, ss2_all = summary_table(f1, alpha=0.05)
        # ss2: Obs, Dep Var Population, Predicted Value, Std Error Mean Predict,
        # Mean ci 95% low, Mean ci 95% upp, Predict ci 95% low, Predict ci 95% upp,
        # Residual, Std Error Residual, Student Residual, Cook's D

        fittedvalues = data[:,2]
        predict_mean_se = data[:,3]
        predict_mean_ci_low, predict_mean_ci_upp = data[:,4:6].T
        predict_ci_low, predict_ci_upp = data[:,6:8].T

        fittedvalues_all = data_all[:,2]
        predict_mean_se_all = data_all[:,3]
        predict_mean_ci_low_all, predict_mean_ci_upp_all = data_all[:,4:6].T
        predict_ci_low_all, predict_ci_upp_all = data_all[:,6:8].T

        for j, kval in enumerate(KindList):
            if kval == 'macro':

                macCiH.append(predict_mean_ci_upp[j])
                macCiL.append(predict_mean_ci_low[j])
                MacPIx.append(Nlist[j])
                MacFitted.append(f.fittedvalues[j])

            elif kval == 'micro':

                micCiH.append(predict_mean_ci_upp[j])
                micCiL.append(predict_mean_ci_low[j])
                MicPIx.append(Nlist[j])
                MicFitted.append(f.fittedvalues[j])

        MicPIx, MicFitted, micCiH, micCiL = zip(*sorted(zip(MicPIx, MicFitted, micCiH, micCiL)))
        MacPIx, MacFitted, macCiH, macCiL = zip(*sorted(zip(MacPIx, MacFitted, macCiH, macCiL)))

        Nlist_all = list(Nlist)
        Nlist_all, predict_mean_ci_low_all, predict_mean_ci_upp_all = zip(*sorted(zip(Nlist_all, predict_mean_ci_low_all, predict_mean_ci_upp_all)))

        num = min(len(MacListX), len(MicListX))
        for i in range(num):
            plt.scatter(MacListX[i], MacListY[i], color = 'LightCoral', alpha= 1 , s = 4, linewidths=0.5, edgecolor='Crimson')
            plt.scatter(MicListX[i], MicListY[i], color = 'SkyBlue', alpha= 1 , s = 4, linewidths=0.5, edgecolor='Steelblue')

        plt.fill_between(MacPIx, macCiL, macCiH, color='r', lw=0.0, alpha=0.3)
        plt.fill_between(MicPIx, micCiL, micCiH, color='b', lw=0.0, alpha=0.3)
        plt.fill_between(Nlist_all, predict_mean_ci_low_all, predict_mean_ci_upp_all, color='0.3', lw=0.0, alpha=0.6)

        MicInt = round(np.mean(MicIntList), 2)
        MicCoef = round(np.mean(MicCoefList), 2)
        MacInt = round(np.mean(MacIntList), 2)
        MacCoef = round(np.mean(MacCoefList), 2)
        R2 = round(np.mean(R2List), 3)

        if index == 0:
            plt.ylim(-0.1, 2.0)
            plt.xlim(1, 6)
            plt.text(1.35, 1.7, r'$micro$'+ ' = '+str(round(MicInt,2))+'+'+str(round(MicCoef, 2))+'*'+r'$N$', fontsize=fs-1, color='Steelblue')
            plt.text(1.35, 1.5, r'$macro$'+ ' = '+str(round(MacInt,2))+'+'+str(round(MacCoef,2))+'*'+r'$N$', fontsize=fs-1, color='Crimson')
            plt.text(1.35, 1.2,  r'$R^2$' + '=' +str(R2), fontsize=fs-1, color='k')

            plt.scatter([0],[-1], color = 'SkyBlue', alpha = 1, s=10, linewidths=0.9, edgecolor='Steelblue', label= 'microbes (n='+str(len(MicListY))+')')
            plt.scatter([0],[-1], color = 'LightCoral',alpha= 1, s=10, linewidths=0.9, edgecolor='Crimson', label= 'macrobes (n='+str(len(MacListY))+')')
            plt.legend(bbox_to_anchor=(-0.04, 1.1, 2.5, .2), loc=10, ncol=2, mode="expand",prop={'size':fs+2})

        elif index == 1:
            plt.ylim(0, 6)
            plt.xlim(1, 6)

            plt.text(1.35, 5.1, r'$micro$'+ ' = '+str(round(MicInt,2))+'+'+str(round(MicCoef,2))+'*'+r'$N$', fontsize=fs-1, color='Steelblue')
            plt.text(1.35, 4.5, r'$macro$'+ ' = '+str(round(MacInt,2))+'+'+str(round(MacCoef,2))+'*'+r'$N$', fontsize=fs-1, color='Crimson')
            plt.text(1.35, 3.75,  r'$R^2$' + '=' +str(R2), fontsize=fs-1, color='k')

        elif index == 2:
            plt.ylim(-2.5, 0.0)
            plt.xlim(1, 6)

            plt.text(1.35, -1.8, r'$micro$'+ ' = '+str(round(MicInt,2))+' '+str(round(MicCoef,2))+'*'+r'$N$', fontsize=fs-1, color='Steelblue')
            plt.text(1.35, -2.0, r'$macro$'+ ' = '+str(round(MacInt,2))+' '+str(round(MacCoef,2))+'*'+r'$N$', fontsize=fs-1, color='Crimson')
            plt.text(1.35, -2.3,  r'$R^2$' + '=' +str(R2), fontsize=fs-1, color='k')

        elif index == 3:
            plt.ylim(0.9, 4.5)
            plt.xlim(1, 6)

            plt.text(1.35, 3.9, r'$micro$'+ ' = '+str(round(MicInt,2))+'+'+str(round(MicCoef,2))+'*'+r'$N$', fontsize=fs-1, color='Steelblue')
            plt.text(1.35, 3.5, r'$macro$'+ ' = '+str(round(MacInt,2))+'+'+str(round(MacCoef,2))+'*'+r'$N$', fontsize=fs-1, color='Crimson')
            plt.text(1.35, 3.0,  r'$R^2$' + '=' +str(R2), fontsize=fs-1, color='k')

        plt.xlabel('Total abundance, ' + r'$log_{10}$', fontsize=fs)
        plt.ylabel(metric, fontsize=fs)
        plt.tick_params(axis='both', which='major', labelsize=fs-3)

    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    #plt.savefig(mydir+'/figs/Fig1/Locey_Lennon_2015_Fig1-OpenReference_NoMicrobeSingletons.png', dpi=600, bbox_inches = "tight")
    #plt.savefig(mydir+'/figs/Fig1/Locey_Lennon_2015_Fig1-ClosedReference_NoMicrobeSingletons.png', dpi=600, bbox_inches = "tight")
    plt.savefig(mydir+'/figs/Fig1/Locey_Lennon_2015_Fig1-ClosedReference.png', dpi=600, bbox_inches = "tight")
    #plt.savefig(mydir+'/figs/Fig1/Locey_Lennon_2015_Fig1-OpenReference.png', dpi=600, bbox_inches = "tight")

    #plt.show()
    plt.close()

    return


Fig1()
