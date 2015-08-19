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


    metlist = []
    metrics = ['Rarity', 'Dominance', 'Evenness','Richness']
    for metric in metrics:

        fig = plt.figure()
        fs = 10 # font size used across figures

        IntList, CoefList, R2List, metlist, KindList, AvgNList, AvgSList = [[], [], [], [], [], [], []]
        Nlist, Slist, Evarlist, ESimplist, klist, radDATA, BPlist, NmaxList, rareSkews, StdList = [[], [], [], [], [], [], [], [], [], []]

        for dataset in datasets:
            f = list()

            its = 1
            for n in range(its):
                radDATA = []

                name, kind, numlines = dataset
                lines = []
                lines = np.random.choice(range(1, numlines+1), numlines, replace=False)

                #path = mydir+'data/'+kind+'/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt'
                path = mydir+'data/'+kind+'/'+name+'/'+name+'-SADMetricData.txt'

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

                    BPlist.append(float(BP))
                    NmaxList.append(float(np.log10(float(Nmax))))

                    # log-modulo transformation of skewnness
                    lms = np.log10(np.abs(float(skew)) + 1)
                    if skew < 0: lms = lms * -1
                    rareSkews.append(float(lms))

                if metric == 'Rarity': metlist = rareSkews
                elif metric == 'Dominance': metlist = NmaxList
                elif metric == 'Evenness': metlist = ESimplist
                elif metric == 'Richness': metlist = Slist

                # Simple regression
                d = pd.DataFrame({'N': list(Nlist)})
                d['y'] = list(metlist)
                f = smf.ols('y ~ N', d).fit()

                IntList.append(f.params[0])
                CoefList.append(f.params[1])
                R2List.append(f.rsquared)
                KindList.append(kind)
                AvgNList.append(np.mean(Nlist))
                AvgSList.append(np.mean(Slist))

        for i in range(4):

            fig.add_subplot(2, 2, i+1)

            if i < 2:
                xlabel = 'Mean total abundance, ' + r'$log_{10}$'
                Xlist = list(AvgNList)
            else:
                xlabel = 'Mean richness, ' + r'$log_{10}$'
                Xlist = list(AvgSList)

            if i == 0 or i == 2:
                Ylist = list(IntList)
                ylabel = 'OLS intercept'
            else:
                Ylist = list(CoefList)
                ylabel = 'OLS coefficient'

            # correlation
            r, r_pval = stats.pearsonr(Xlist, Ylist)
            rho, rho_pval = stats.spearmanr(Xlist, Ylist)

            micxlist = []
            micylist = []
            macxlist = []
            macylist = []

            for k, val in enumerate(KindList):
                if val == 'micro':
                    micxlist.append(Xlist[k])
                    micylist.append(Ylist[k])

                elif val == 'macro':
                    macxlist.append(Xlist[k])
                    macylist.append(Ylist[k])

            plt.scatter(micxlist, micylist, color = 'SkyBlue', alpha= 1,
                        s = 20, linewidths=0.9, edgecolor='Steelblue')

            plt.scatter(macxlist, macylist, color = 'LightCoral', alpha= 1,
                        s = 20, linewidths=0.9, edgecolor='Crimson')

            #elif i == 1 or i == 3: plt.ylim(min(Ylist), max(Ylist))
            #elif i == 0 or i == 2: plt.ylim(min(Ylist), max(Ylist))

            txt = r'$r$'+' = '+str(round(r,2))+', '+'$p$'+' = '+str(round(r_pval,3))
            txt = txt +'\n'+r'$rho$'+' = '+str(round(rho,2))+', '+'$p$'+' = '+str(round(rho_pval,3))

            plt.scatter([0],[min(Ylist)], alpha = 0, s=0, label=txt)
            leg = plt.legend(loc=1, prop={'size':fs-1}, numpoints=1)
            leg.draw_frame(False)

            plt.xlim(min(Xlist), max(Xlist)*1.5)
            plt.xlabel(xlabel, fontsize=fs)
            plt.ylabel(ylabel, fontsize=fs)
            plt.tick_params(axis='both', which='major', labelsize=fs-3)

        print metric
        plt.subplots_adjust(wspace=0.4, hspace=0.4)
        #plt.savefig(mydir+'/figs/appendix/'+metric+'.png', dpi=600, bbox_inches = "tight")
        #sys.exit()
        plt.show()

    plt.close()

    return


Fig1()
