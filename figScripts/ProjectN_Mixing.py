from __future__ import division
import  matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
import random
import scipy as sc

import sys

sys.path.append("/Users/lisalocey/Desktop/RareBio/")
import feasible_functions as ff

radDATA = open('/Users/lisalocey/Desktop/RareBio/RADdataEMP.txt','r')
#radDATA = open('/Users/lisalocey/Desktop/RareBio/RADdataECON.txt','r')


fs = 10
color = str()

Nlist, Slist, Evarlist, ESimplist, ENeelist, EHeiplist, EQlist = [], [], [], [], [], [], []
klist, Shanlist, BPlist, SimpDomlist, SinglesList, tenlist, onelist = [], [], [], [], [], [], []

ptonelist, ptzonelist, BP1list, BP2list, BP3list, BP4list, BP5list = [], [], [], [], [], [], []
NmaxList, rareOnesList, rareRelList, rarePairList, rareSumOnesList = [], [], [], [], []

ct = 0
klist = []
for data in radDATA:
    
    data_list = data.split()
    #print 'datalist', data_list
    name, kind, N, S, Evar, ESimp, ENee, EHeip, EQ, EPielou, BP, SimpDom, perc_ones, s_80, s_ten, s_one, s_ptone, s_ptzone, BP1, BP2, BP3, BP4, BP5, Nmax, Nperc90, Nperc70, Nperc50, Nperc30, Nperc10, rareRel, rareOnes, rarePair, rareSumOnes = data_list
    
    N = float(N)
    S = float(S)
    
    Nlist.append(N)
    Slist.append(S)
    
    Evarlist.append(float(Evar))
    ESimplist.append(float(ESimp))
    ENeelist.append(float(ENee))
    Shanlist.append(float(EPielou))
    
    BPlist.append(float(BP))
    SimpDomlist.append(float(SimpDom))
    
    EHeiplist.append(float(EHeip))
    EQlist.append(float(EQ))
    tenlist.append(float(s_ten))
    onelist.append(float(s_one))
    ptonelist.append(float(s_ptone))
    ptzonelist.append(float(s_ptzone))
    SinglesList.append(float(perc_ones))
    
    BP1list.append(float(BP1))
    BP2list.append(float(BP2))
    BP3list.append(float(BP3))
    BP4list.append(float(BP4))
    BP5list.append(float(BP5))
    NmaxList.append(float(Nmax))
    
    rareRelList.append(float(rareRel))
    rareOnesList.append(float(rareOnes))
    rarePairList.append(float(rarePair))
    rareSumOnesList.append(float(rareSumOnes))
    
    if kind == 'micro':
        color = 'm'
        
    elif kind == 'macro':
        color = 'c'
    
    else: color = 'Steelblue'
    
    klist.append(color)
  
    
    ct+=1
    #if ct >= 1000: break

    
for i, val in enumerate(Nlist):
    if val == NmaxList[i]:
        print 'most abundant = N'
    if val < NmaxList[i] + 9:
        print 'S cannot be >= 10'
            
            

    
 

def scatterNmax():
    
    metrics = [['Most abundant species', NmaxList]]
    
    for i in metrics:
    
        fig = plt.figure()
        ax = fig.add_subplot(2, 2, 1)

        metric = i[0]
        metlist = i[1]
        
        nmacs = 0
        nmics = 0
        
        for j, k in enumerate(klist):
            
            plt.scatter(Nlist[j], metlist[j], color = k, alpha = 0.6, s = 2,
            linewidths = 0.1)
            
            if k == 'm': nmics += 1
            elif k == 'c': nmacs += 1
            
        plt.scatter([0],[-1], color = 'm', alpha = 0.8, s=10,
                     label= 'microbial (n='+str(nmics)+')')
        
        plt.scatter([0],[-1], color = 'c',alpha=0.8, s=10,
                     label= 'macrobial (n='+str(nmacs)+')')
        
        plt.legend(bbox_to_anchor=(-0.03, 1.1, 2.47, .2), loc=10, ncol=2,
                               mode="expand",prop={'size':fs+2})
        
        
        xmax = max(Nlist)
        plt.plot([1,xmax],[1, xmax], 'k--', lw=1)
        
        X = list(np.log(Nlist))
        Y = list(np.log(metlist))
        slope, intercept, rval, pval, stderr = sc.stats.linregress(X,Y)
        print 'r-squared and slope for scatterNmax:',rval**2, slope
        
        plt.xscale('log')
        plt.yscale('log')
        
        plt.xlabel('Total abundance', fontsize=fs)
        plt.ylabel(metric, fontsize=fs)
        plt.tick_params(axis='both', which='major', labelsize=fs-3)
        
        ax = fig.add_subplot(2, 2, 2)
    
        for j, k in enumerate(klist):
            
            plt.scatter(Slist[j], metlist[j], color = k, alpha = 0.6, s = 2,
            linewidths = 0.1)
            
        
        X = list(np.log(Slist))
        Y = list(np.log(metlist))
        slope, intercept, rval, pval, stderr = sc.stats.linregress(X,Y)
        print 'r-squared and slope for scatterNmax:',rval**2, slope
        
        #z = np.polyfit(X,Y,1)
        #p = np.poly1d(z)
        #xp = np.linspace(min(X), max(X), 1000)
        #plt.plot(xp,p(xp),'-',c='k',lw=1)

        plt.xlim(min(Slist), max(Slist))
        plt.ylim(min(metlist), max(metlist))
        
        plt.xlabel('Richness', fontsize=fs)
        plt.ylabel(metric, fontsize=fs)
        plt.tick_params(axis='both', which='major', labelsize=fs-3)
        
        plt.xscale('log')
        plt.yscale('log')
            
    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/ProjectN/ScatterMax.png',
    dpi=600)
    plt.close()
        
    return



    
#scatterMetrics()
#scatterMetricsRare()
#kdensMetrics()
#scatterNmax()
#scatterNmax2()
#MoneyFig1()
#MoneyFig2()