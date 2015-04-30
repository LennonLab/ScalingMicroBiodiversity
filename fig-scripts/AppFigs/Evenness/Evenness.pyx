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
    BadNames = ['.DS_Store', 'EMPclosed', 'EMPopen', 'BCI', 'AGSOIL', 'SLUDGE', 'FECES']

    for name in os.listdir(mydir2 +'data/micro'):
        if name in BadNames: continue

        path = mydir2+'data/micro/'+name+'/'+name+'-SADMetricData.txt' # _NoMicrobe1s
        num_lines = sum(1 for line in open(path))
        datasets.append([name, 'micro', num_lines])

    for name in os.listdir(mydir2 +'data/macro'):
        if name in BadNames: continue

        path = mydir2+'data/macro/'+name+'/'+name+'-SADMetricData.txt' # _NoMicrobe1s
        num_lines = sum(1 for line in open(path))
        datasets.append([name, 'macro', num_lines])


    metrics = ['Evar', 'ESimp, '+r'$log_{10}$', 'Heip\'s', 'O']

    fig = plt.figure()
    for index, i in enumerate(metrics):

        metric = i
        fig.add_subplot(2, 2, index+1)
        fs = 10 # font size used across figures

        MicIntList, MicCoefList, MacIntList, MacCoefList = [[], [], [], []]
        Nlist, Slist, ESimplist, klist, radDATA, BPlist, NmaxList, rareSkews, KindList = [[], [], [], [], [], [], [], [], []]
        PrestonAList, EvarList, EQList, EHeipList, OList = [[],[],[],[],[]]
        SimpDomList, McNList, LogSkewList, POnesList = [[],[],[],[]]

        its = 100
        for n in range(its):
            print n, metric

            Nlist, Slist, ESimplist, klist, radDATA, BPlist, NmaxList, rareSkews, KindList = [[], [], [], [], [], [], [], [], []]
            PrestonAList, EvarList, EQList, EHeipList, OList = [[],[],[],[],[]]
            SimpDomList, McNList, LogSkewList, POnesList = [[],[],[],[]]

            numMac = 0
            numMic = 0
            radDATA = []

            for dataset in datasets:

                name, kind, numlines = dataset
                lines = []

                if numlines > 40: lines = random.sample(range(1, numlines+1), 40)
                else: lines = random.sample(range(1, numlines+1), 40)

                path = mydir2+'data/'+kind+'/'+name+'/'+name+'-SADMetricData.txt'

                for line in lines:
                    data = linecache.getline(path, line)
                    radDATA.append(data)

            for data in radDATA:

                data = data.split()
                name, kind, N, S, Evar, Camargo, ESimp, EQ, O, Camargo, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, p_ones, p_zpt1, preston = data

                KindList.append(kind)
                N = float(N)
                S = float(S)

                if S < 2 or N < 10: continue # Min species richness

                Nlist.append(float(np.log10(N)))
                Slist.append(float(np.log10(S)))

                # Evenness
                ESimplist.append(float(np.log10(float(ESimp))))
                EvarList.append(float(np.log10(float(Evar))))
                EHeipList.append(float(np.log10(float(EHeip))))
                OList.append(float(np.log10(float(O))))

                PrestonAList.append(float(preston))

            MacListX = []
            MacListY = []
            MicListX = []
            MicListY = []
            metlist = []

            if index == 0: metlist = list(EvarList)
            elif index == 1: metlist = list(ESimplist)
            elif index == 2: metlist = list(EHeipList)
            elif index == 3: metlist = list(OList)

            for j, k in enumerate(KindList):
                if k == 'micro':
                    MicListX.append(Nlist[j])
                    MicListY.append(metlist[j])

                elif k == 'macro':
                    MacListX.append(Nlist[j])
                    MacListY.append(metlist[j])

            # Multiple regression
            d = pd.DataFrame({'N': list(Nlist)})
            d['y'] = list(metlist)
            d['Kind'] = list(KindList)

            f = smf.ols('y ~ N * Kind', d).fit()

            MacIntList.append(f.params[0])
            MacCoefList.append(f.params[2])
            MicIntList.append(f.params[1] + f.params[0])
            MicCoefList.append(f.params[3] + f.params[2])

            r2 = f.rsquared
            pvals = f.pvalues

        MacPIx, MacFitted, MicPIx, MicFitted = [[],[],[],[]]
        macCiH, macCiL, micCiH, micCiL = [[],[],[],[]]

        lm = smf.ols('y ~ N * Kind', d).fit()
        print metric, '\n', lm.summary()
        f1 = smf.ols('y ~ N', d).fit()
        print metric, '\n', f1.summary()

        st, data, ss2 = summary_table(lm, alpha=0.05)
        # ss2: Obs, Dep Var Population, Predicted Value, Std Error Mean Predict,
        # Mean ci 95% low, Mean ci 95% upp, Predict ci 95% low, Predict ci 95% upp,
        # Residual, Std Error Residual, Student Residual, Cook's D

        fittedvalues = data[:,2]
        predict_mean_se = data[:,3]
        predict_mean_ci_low, predict_mean_ci_upp = data[:,4:6].T
        predict_ci_low, predict_ci_upp = data[:,6:8].T

        MacInt = lm.params[0]
        MacCoef = lm.params[2]
        MicInt = lm.params[1] + MacInt
        MicCoef = lm.params[3] + MacCoef


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

        for i in range(len(MicListX)):
            plt.scatter(MacListX[i], MacListY[i], color = 'LightCoral', alpha= 1 , s = 4, linewidths=0.5, edgecolor='Crimson')
            plt.scatter(MicListX[i], MicListY[i], color = 'SkyBlue', alpha= 1 , s = 4, linewidths=0.5, edgecolor='Steelblue')

        plt.fill_between(MacPIx, macCiL, macCiH, color='r', lw=0.0, alpha=0.3)
        plt.fill_between(MicPIx, micCiL, micCiH, color='b', lw=0.0, alpha=0.3)

        MicInt = round(np.mean(MicIntList), 2)
        MicCoef = round(np.mean(MicCoefList), 2)
        MacInt = round(np.mean(MacIntList), 2)
        MacCoef = round(np.mean(MacCoefList), 2)

        if index == 0:
            plt.ylim(-1.4, 0.1)
            #plt.xlim(1, 7)
            plt.text(.2, -1.1, r'$y_{micro}$'+ ' = '+str(round(MacInt,2))+'+'+str(round(MacCoef,2))+'*'+r'$N$', fontsize=fs-1, color='Steelblue')
            plt.text(.2, -1.22, r'$y_{macro}$'+ ' = '+str(round(MacInt,2))+'+'+str(round(MacCoef,2))+'*'+r'$N$', fontsize=fs-1, color='Crimson')
            plt.text(.2, -1.33,  r'$R^2$' + '=' +str(round(r2,3)), fontsize=fs-1, color='k')


        if index == 1:
            #plt.ylim(0, 100)
            #plt.xlim(1, 8)
            plt.text(.2, -1.8, r'$y_{micro}$'+ ' = '+str(round(MacInt,2))+'+'+str(round(MacCoef,2))+'*'+r'$N$', fontsize=fs-1, color='Steelblue')
            plt.text(.2, -2.1, r'$y_{macro}$'+ ' = '+str(round(MacInt,2))+'+'+str(round(MacCoef,2))+'*'+r'$N$', fontsize=fs-1, color='Crimson')
            plt.text(.2, -2.4,  r'$R^2$' + '=' +str(round(r2,3)), fontsize=fs-1, color='k')

        if index == 2:
            #plt.ylim(0, 1)
            #plt.xlim(1, 7)
            plt.text(.2, -1.3, r'$y_{micro}$'+ ' = '+str(round(MacInt,2))+'+'+str(round(MacCoef,2))+'*'+r'$N$', fontsize=fs-1, color='Steelblue')
            plt.text(.2, -1.5, r'$y_{macro}$'+ ' = '+str(round(MacInt,2))+'+'+str(round(MacCoef,2))+'*'+r'$N$', fontsize=fs-1, color='Crimson')
            plt.text(.2, -1.8,  r'$R^2$' + '=' +str(round(r2,3)), fontsize=fs-1, color='k')

        if index == 3:
            plt.ylim(-1.4, 0)
            #plt.xlim(1, 7)
            plt.text(.2, -1.0, r'$y_{micro}$'+ ' = '+str(round(MacInt,2))+'+'+str(round(MacCoef,2))+'*'+r'$N$', fontsize=fs-1, color='Steelblue')
            plt.text(.2, -1.12, r'$y_{macro}$'+ ' = '+str(round(MacInt,2))+'+'+str(round(MacCoef,2))+'*'+r'$N$', fontsize=fs-1, color='Crimson')
            plt.text(.2, -1.24,  r'$R^2$' + '=' +str(round(r2,3)), fontsize=fs-1, color='k')


        plt.xlabel('Total abundance, ' + r'$log_{10}$', fontsize=fs-2)
        plt.ylabel(metric, fontsize=fs-2)
        plt.tick_params(axis='both', which='major', labelsize=fs-3)

    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    plt.savefig(mydir+'/figs/appendix/AppendixEvennessFig_NoEMP.png', dpi=600, bbox_inches = "tight")
    plt.close()

    return


Fig1()
