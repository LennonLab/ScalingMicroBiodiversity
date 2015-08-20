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
    #BadNames = ['.DS_Store', 'BCI', 'AGSOIL', 'SLUDGE', 'EMPopen', 'BOVINE', 'FECES', 'MGRASTopen', 'MGRAST', 'NABC', 'FUNGI']

    for name in os.listdir(mydir +'data/micro'):
        if name in GoodNames: pass
        else: continue

        #if name in BadNames: continue
        #else: pass

        #path = mydir+'data/micro/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt'
        path = mydir+'data/micro/'+name+'/'+name+'-SADMetricData.txt'

        num_lines = sum(1 for line in open(path))
        datasets.append([name, 'micro', num_lines])
        print name, num_lines

    for name in os.listdir(mydir +'data/macro'):
        if name in GoodNames: pass
        else: continue

        #if name in BadNames: continue
        #else: pass

        #path = mydir+'data/macro/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt'
        path = mydir+'data/macro/'+name+'/'+name+'-SADMetricData.txt'

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
                if numlines > 1000:
                    lines = np.random.choice(range(1, numlines+1), 1000, replace=False)
                else:
                    lines = np.random.choice(range(1, numlines+1), numlines, replace=False)

                #path = mydir+'data/'+kind+'/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt'
                path = mydir+'data/'+kind+'/'+name+'/'+name+'-SADMetricData.txt'

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

                if f.params[0] > 1.0:
                    print name

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
                ylabel = 'Intercept'
            else:
                Ylist = list(CoefList)
                ylabel = 'Scaling exponent'

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


            txt = r'$r$'+' = '+str(round(r,2))+', '+'$p$'+' = '+str(round(r_pval,3))
            #txt = txt +'\n'+r'$rho$'+' = '+str(round(rho,2))+', '+'$p$'+' = '+str(round(rho_pval,3))

            plt.scatter(np.mean(Xlist),np.mean(Ylist), alpha = 0, s=0, label=txt)

            plt.legend(bbox_to_anchor=(-0.01, 1.01, 1.02, .2), loc=10, ncol=1,
                                mode="expand",prop={'size':fs}, numpoints=1)

            titletext = "The relationship of "+metric+" to sample abundance should not be\n"
            titletext += "influenced by sample abundance. Blue dots are microbial datasets.\n"
            titletext += "Red dots are datasets of trees, birds, and mammals."

            if i == 0:
                fig.suptitle(titletext, y=1.07)

            #plt.ylim(0.95*min(Ylist), 1.05*max(Ylist))
            #plt.xlim(0.95*min(Xlist), 1.05*max(Xlist))
            plt.xlabel(xlabel, fontsize=fs)
            plt.ylabel(ylabel, fontsize=fs)
            plt.tick_params(axis='both', which='major', labelsize=fs-3)

        print metric
        plt.subplots_adjust(wspace=0.4, hspace=0.5)
        #plt.savefig(mydir+'/figs/appendix/OLSmodel_Checks_PerDataset/'+metric+'_OpenRef.png', dpi=600, bbox_inches = "tight")
        #plt.savefig(mydir+'/figs/appendix/OLSmodel_Checks_PerDataset/'+metric+'_OpenRef_NoMicrobe1s.png', dpi=600, bbox_inches = "tight")
        plt.savefig(mydir+'/figs/appendix/OLSmodel_Checks_PerDataset/'+metric+'_ClosedRef.png', dpi=600, bbox_inches = "tight")
        #plt.savefig(mydir+'/figs/appendix/OLSmodel_Checks_PerDataset/'+metric+'_ClosedRef_NoMicrobe1s.png', dpi=600, bbox_inches = "tight")
        #plt.show()

    #plt.close()

    return


Fig1()
