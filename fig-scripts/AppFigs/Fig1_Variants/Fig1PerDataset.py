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


def Fig1():

    OUT = open(mydir + 'output/PerDataset.txt','w+')

    """ This code generates a 4 plot figure of diversity properties (rarity,
        dominance, evenness, richness) versus total abundance, for each dataset.
        This code also generates a .txt file of results for the regression
        analyses. """

    datasets = []
    GoodNames = ['MGRAST', 'HMP', 'EMPclosed', 'BBS', 'CBC', 'MCDB', 'GENTRY', 'FIA']

    for name in os.listdir(mydir +'data/micro'):
        if name in GoodNames: pass
        else: continue

        path = mydir+'data/micro/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt'
        #path = mydir+'data/micro/'+name+'/'+name+'-SADMetricData.txt'

        num_lines = sum(1 for line in open(path))
        datasets.append([name, 'micro', num_lines])
        print name, num_lines

    for name in os.listdir(mydir +'data/macro'):
        if name in GoodNames: pass
        else: continue

        path = mydir+'data/macro/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt'
        #path = mydir+'data/macro/'+name+'/'+name+'-SADMetricData.txt'

        num_lines = sum(1 for line in open(path))
        datasets.append([name, 'macro', num_lines])
        print name, num_lines


    metrics = ['Rarity, '+r'$log_{10}$',
            'Dominance, '+r'$log_{10}$',
            'Evenness, ' +r'$log_{10}$',
            'Richness, ' +r'$log_{10}$']

    OUT = open(mydir + 'output/SummaryPerDataset_NoMicrobe1s.txt','w+')
    #OUT = open(mydir + 'output/SummaryPerDataset.txt','w+')

    for dataset in datasets:

        fig = plt.figure()
        for index, i in enumerate(metrics):

            metric = i
            fig.add_subplot(2, 2, index+1)
            fs = 10 # font size used across figures

            IntList, CoefList, R2List, metlist = [[], [], [], []]
            Nlist, Slist, Evarlist, ESimplist, klist, radDATA, BPlist, NmaxList, rareSkews, KindList, StdList = [[], [], [], [], [], [], [], [], [], [], []]
            #name, kind, N, S, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]

            its = 1000
            f = list()
            for n in range(its):

                #name, kind, N, S, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
                Nlist, Slist, Evarlist, ESimplist, klist, radDATA, BPlist, NmaxList, rareSkews, KindList, StdList = [[], [], [], [], [], [], [], [], [], [], []]

                numMac = 0
                numMic = 0
                radDATA = []

                name, kind, numlines = dataset
                lines = []
                lines = np.random.choice(range(1, numlines+1), 100, replace=True)

                path = mydir+'data/'+kind+'/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt'
                #path = mydir+'data/'+kind+'/'+name+'/'+name+'-SADMetricData.txt'

                for line in lines:
                    data = linecache.getline(path, line)
                    radDATA.append(data)

                for data in radDATA:

                    data = data.split()
                    name, kind, N, S, Var, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = data

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

                if index == 0: metlist = list(rareSkews)
                elif index == 1: metlist = list(NmaxList)
                elif index == 2: metlist = list(ESimplist)
                elif index == 3: metlist = list(Slist)

                # Simple regression
                d = pd.DataFrame({'N': list(Nlist)})
                d['y'] = list(metlist)
                f = smf.ols('y ~ N', d).fit()

                IntList.append(f.params[0])
                CoefList.append(f.params[1])
                R2List.append(f.rsquared)


            PIx = list(Nlist)
            st, data, ss2 = summary_table(f, alpha=0.05)
            # ss2: Obs, Dep Var Population, Predicted Value, Std Error Mean Predict,
            # Mean ci 95% low, Mean ci 95% upp, Predict ci 95% low, Predict ci 95% upp,
            # Residual, Std Error Residual, Student Residual, Cook's D

            Fitted = data[:,2]
            predict_mean_se = data[:,3]
            CiL, CiH = data[:,4:6].T
            PiL, PiH = data[:,6:8].T

            PIx, Fitted, CiH, CiL = zip(*sorted(zip(PIx, Fitted, CiH, CiL)))

            plt.scatter(Nlist, metlist, color = 'SkyBlue', alpha= 1 , s = 4, linewidths=0.5, edgecolor='Steelblue')
            plt.fill_between(PIx, CiL, CiH, color='b', lw=0.0, alpha=0.3)

            Int = round(np.mean(IntList), 2)
            Coef = round(np.mean(CoefList), 2)
            R2 = round(np.mean(R2List), 3)

            print dataset, metric, Int, Coef, R2

            x = min(Nlist)
            y = 1.1*max(metlist)

            plt.scatter([0],[-1], color = 'SkyBlue', alpha = 1, s=10, linewidths=0.9, edgecolor='Steelblue', label= metric+' = '+str(round(Int,2))+'+'+str(round(Coef, 2))+'*'+r'$N$'+'\n'+r'$R^2$' + '=' +str(R2) +' (n='+str(len(PIx))+')')

            if index == 2:
                leg = plt.legend(loc=3,prop={'size':fs-1})
                leg.draw_frame(False)

            else:
                leg = plt.legend(loc=2,prop={'size':fs-1})
                leg.draw_frame(False)

            plt.ylim(min(metlist), max(metlist)*1.1)
            plt.xlim(min(Nlist), max(Nlist))

            plt.xlabel('Total abundance, ' + r'$log_{10}$', fontsize=fs)
            plt.ylabel(metric, fontsize=fs)
            plt.tick_params(axis='both', which='major', labelsize=fs-3)

            metrix = ['rarity', 'dominance', 'evenness', 'richness']
            print>>OUT, name, kind, metrix[index], np.mean(PIx), np.mean(Slist), Int, Coef

        #plt.subplots_adjust(wspace=0.4, hspace=0.4)
        plt.savefig(mydir+'/figs/appendix/Fig1/PerDataset/Locey_Lennon_2015_'+name+'_NoMicrobeSingletons.png', dpi=600, bbox_inches = "tight")
        #plt.show()
        print name

        plt.close()
    OUT.close()
    return


Fig1()
