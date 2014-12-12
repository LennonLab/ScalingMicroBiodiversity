from __future__ import division
import  matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
import random
import scipy as sc
import os
import sys

import pandas
from pandas.tools import plotting
from scipy import stats
from statsmodels.formula.api import ols
from numpy.random import randn

from math import log10


mydir = os.path.expanduser("~/Desktop/Repos/rare-bio/")

sys.path.append(mydir + "tools/feasible_functions")
import feasible_functions as ff

sys.path.append(mydir + "tools/partitions")
#import partitions


OrC = 'open'


def get_kdens(summands):
    """ Finds the kernel density function across a sample of parts
    of partitions for a given total (N) and number of parts (S) """
    
    density = gaussian_kde(summands)
    n = 1000 #len(summands)
    xs = np.linspace(float(min(summands)),float(max(summands)),n)
    density.covariance_factor = lambda : .4
    density._compute_covariance()
    D = [xs,density(xs)]
    return D


fs = 10 # font size used across figures
color = str()



def Fig1():

    Nlist, Slist, Evarlist, ESimplist, ENeelist, EHeiplist, EQlist = [[], [], [], [], [], [], []]
    klist, Shanlist, BPlist, SimpDomlist, SinglesList, tenlist, onelist = [[], [], [], [], [], [], []]
    NmaxList, rareOnesList, rareRelList, rarePairList, rareSkews, KindList = [[], [], [], [], [], []]
    NSlist = []
    
    ct = 0
    klist = []
    
    radDATA = ff.radDATA()
    
    for data_list in radDATA:
        
        name, kind, N, S, Evar, ESimp, ENee, EHeip, EQ, EPielou, BP, SimpDom, rareRel, rareOnes, skew = data_list
        
        N = float(N)
        S = float(S)
        
        if S < 4: continue
        
        N = float(N)
        S = float(S)
        
        Nlist.append(float(np.log(N)))
        Slist.append(float(np.log(S)))
        NSlist.append(float(np.log(N/S)))
        
        Evarlist.append(float(np.log(float(Evar))))
        ESimplist.append(float(np.log(float(ESimp))))
        KindList.append(kind)
        
        BPlist.append(float(BP))
        NmaxList.append(float(np.log(float(BP)*float(N))))
        EHeiplist.append(float(EHeip))
        
        rareRelList.append(float(rareOnes/S))
        rareOnesList.append(float(rareOnes))
        rareSkews.append(float(skew))
        
        if kind == 'micro': klist.append('b')
        if kind == 'macro': klist.append('r')
        
        ct+=1
    
        
    
    metrics = [['log(skewnness)', rareSkews], 
              ['log(greatest species abundance)', NmaxList], 
              ['log(species evenness)', Evarlist]]
            


    fig = plt.figure()
    for index, i in enumerate(metrics):
    
        ax = fig.add_subplot(2, 3, index+1)
        ax.set_axis_bgcolor('w')
        
        metric = i[0]
        metlist = i[1]
        
        MacListX = []
        MacListY = []
        
        MicListX = []
        MicListY = []
        
        nmacs = 0
        nmics = 0
        
        for j, k in enumerate(klist):
            
            if k == 'b': 
                nmics += 1
                if index >= 0: MicListX.append(Nlist[j])
                else: MicListX.append(Slist[j])
                
                if metlist[j] == 0: MicListY.append(1)
                else: MicListY.append(metlist[j])
        
                
                
            elif k == 'r':
                nmacs += 1
                if index >= 0: MacListX.append(Nlist[j])
                else: MacListX.append(Slist[j])
                
                if metlist[j] == 0: MacListY.append(1)
                else: MacListY.append(metlist[j])
        
        
        if index == 0:
            #plt.ylim(0, 7)
            plt.xlim(0, 16)
            pass
            
        if index == 1:
            plt.ylim(0, 16)
            plt.xlim(0, 16)
            pass
                
        if index == 2:
            plt.ylim(-3.5, 0)
            plt.xlim(0, 16)
            pass
            
        # scatter plots
        
        indices = range(500)
        random.shuffle(indices)            
        
        headers = ["Category","Skewnness","Evenness","Dominance","N","S","AvgAb"] 
        
        MultipleRegression(KindList, rareSkews, Evarlist, NmaxList, Nlist, Slist, NSlist, headers)
        
        for i in indices:
            plt.scatter(MacListX[i], MacListY[i], color = 'LightCoral', alpha= 1 , s = 4, linewidths=0.5, edgecolor='Crimson')                                                
            plt.scatter(MicListX[i], MicListY[i], color = 'SkyBlue', alpha= 1 , s = 4, linewidths=0.5, edgecolor='Steelblue')
        
        
        Y = np.array(MacListY)
        X = np.array(MacListX)
            
        Macslope, Macintercept, Macrval, Macpval, Macstderr = sc.stats.linregress(X,Y)
        print metric,': r-squared and slope for macrobes:',Macrval**2, Macslope
        
        z = np.polyfit(X,Y,1)
        p = np.poly1d(z)
        xp = np.linspace(min(X)/2, 2*max(X), 1000)
        plt.plot(xp,p(xp),'-',c='r',lw=1)
        
        
        Y = np.array(MicListY)
        X = np.array(MicListX)
            
        Micslope, Micintercept, Micrval, Micpval, Micstderr = sc.stats.linregress(X,Y)
        print metric,': r-squared and slope for microbes:',Micrval**2, Micslope
        
        z = np.polyfit(X,Y,1)
        p = np.poly1d(z)
        xp = np.linspace(min(X)/2, 2*max(X), 1000)
        plt.plot(xp,p(xp),'-',c='b',lw=1)
        
        if index == 0:
            x = 0.5
            y = 6.3
            y2 = 5.8
        elif index == 1:
            x = 0.75
            y = 14.5
            y2 = 13.5
        elif index == 2:
            x = 0.75
            y = -3.1
            y2 = -3.3
        
        
        plt.text(x, y, 'slope ='+str(round(Micslope,2))+', ' + r'$R^2$' + '=' +str(round(Micrval**2,2)), fontsize=fs-2, color='Steelblue')         
        plt.text(x, y2, 'slope ='+str(round(Macslope,2))+', ' + r'$R^2$' + '=' +str(round(Macrval**2,2)), fontsize=fs-2, color='Crimson')         
        
        
        if index == 0:
            plt.scatter([0],[-1], color = 'SkyBlue', alpha = 1, s=10, linewidths=0.9, edgecolor='Steelblue',
                        label= 'microbes (n=500)')
            
            plt.scatter([0],[-1], color = 'LightCoral',alpha= 1, s=10, linewidths=0.9, edgecolor='Crimson',
                        label= 'macrobes (n=500)')
            
            plt.legend(bbox_to_anchor=(-0.03, 1.1, 3.9, .2), loc=10, ncol=2,
                                mode="expand",prop={'size':fs})
        
        
        if index > 0: plt.xlabel('log(abundance of a sample)', fontsize=fs-2)
        else: plt.xlabel('log(number of species)', fontsize=fs-2)
        
        plt.ylabel(metric, fontsize=fs-2)
        plt.tick_params(axis='both', which='major', labelsize=fs-3)
        
    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    plt.savefig(mydir+'/figs/Locey_Lennon_Fig1-'+OrC+'_REF.png', dpi=600, bbox_inches = "tight")
    plt.close()
    
    return




def Fig2():

    Nlist, Slist, Evarlist, ESimplist, ENeelist, EHeiplist, EQlist = [[], [], [],
                                                                    [], [], [], []]
    klist, Shanlist, BPlist, SimpDomlist, SinglesList, tenlist, onelist = [[], [],
                                                                [], [], [], [], []]
    
    NmaxList, rareOnesList, rareRelList, rareSkews = [[], [], [], []]
    
    ct = 0
    klist = []
    
    #SsadDATA = open(mydir+'output/EMP'+OrC+'-SSAD-ResultsTable.txt','r')
    SsadDATA = open(mydir+'output/Macro-SSAD-ResultsTable.txt','r')
    
    for data in SsadDATA:
            
        data_list = data.split()
        N, S, Evar, ESimp, ENee, EHeip, EQ, EPielou, BP, SimpDom, rareRel, rareOnes, skew = data_list
            
        N = float(N)
        S = float(S)
        
        Evarlist.append(float(Evar))    
        Nlist.append(N)
        Slist.append(S)
            
        ESimplist.append(float(ESimp))
        
        BPlist.append(float(BP))
        NmaxList.append(float(BP)*float(N))
        
        EHeiplist.append(float(EHeip))
        rareOnesList.append(float(rareOnes))
        
        rareSkews.append(float(skew))
        
        klist.append('DarkOrange')
        
        ct+=1
        
            
    
    #RadDATA = open(mydir+'output/EMP'+OrC+'-RADdata.txt','r')        
    RadDATA = open(mydir+'output/Macro-RADdata.txt','r')        
    
    
    for data in RadDATA:
            
        data_list = data.split()
        N, S, Evar, ESimp, ENee, EHeip, EQ, EPielou, BP, SimpDom, rareRel, rareOnes, skew = data_list
            
        N = float(N)
        S = float(S)
        
        Evarlist.append(float(Evar))    
        Nlist.append(N)
        Slist.append(S)
            
        ESimplist.append(float(ESimp))
        
        BPlist.append(float(BP))
        NmaxList.append(float(BP)*float(N))
        
        EHeiplist.append(float(EHeip))
        rareOnesList.append(float(rareOnes))
        
        rareSkews.append(float(skew))
        
        klist.append('DarkCyan')
        
        ct+=1
        

        
    metrics = [['log(number of singletons)', rareOnesList], 
              ['log(greatest abundance)', NmaxList], 
              ['log(evenness)', Evarlist]]
            


    fig = plt.figure()
    for index, i in enumerate(metrics):
    
        ax = fig.add_subplot(2, 3, index+1)
        ax.set_axis_bgcolor('w')
        
        metric = i[0]
        metlist = i[1]
        
        RADListX = []
        RADListY = []
        
        SSADListX = []
        SSADListY = []
        
        rads = 0
        ssads = 0
        
        for j, k in enumerate(klist):
            
            if k == 'DarkOrange': 
                ssads += 1
                if index > 0: SSADListX.append(Nlist[j])
                else: SSADListX.append(Slist[j])
                
                if metlist[j] == 0: SSADListY.append(1)
                else: SSADListY.append(metlist[j])
        
                
                
            elif k == 'DarkCyan':
                rads += 1
                if index > 0: RADListX.append(Nlist[j])
                else: RADListX.append(Slist[j])
                
                if metlist[j] == 0: RADListY.append(1)
                else: RADListY.append(metlist[j])
                
        
        if index == 0:
            plt.ylim(0, 9)
            plt.xlim(0, 9)
            pass
            
        if index == 1:
            plt.ylim(0, 16)
            plt.xlim(0, 16)
            pass
                
        if index == 2:
            plt.ylim(-3.5, 0.5)
            plt.xlim(0, 16)
            pass
        
        #print len(SSADListX), len(SSADListY), len(RADListX), len(RADListY)
        #sys.exit()
        
        # scatter plots
        indices = range(1000)
        random.shuffle(indices)            
        
        RADListX = np.log(RADListX)
        RADListY = np.log(RADListY)
        
        SSADListX = np.log(SSADListX)
        SSADListY = np.log(SSADListY)
        
        for i in indices:
            plt.scatter(SSADListX[i], SSADListY[i], color = 'PeachPuff', alpha= 1 , s = 4, linewidths=0.5, edgecolor='DarkOrange')                                                
            plt.scatter(RADListX[i], RADListY[i], color = 'SkyBlue', alpha= 1 , s = 4, linewidths=0.5, edgecolor='DarkCyan')
        
        
        Y = np.array(RADListY)
        X = np.array(RADListX)
            
        RADslope, RADintercept, RADrval, RADpval, RADstderr = sc.stats.linregress(X,Y)
        print metric,': r-squared and slope for RADs:',RADrval**2, RADslope
        
        z = np.polyfit(X,Y,1)
        p = np.poly1d(z)
        xp = np.linspace(min(X)/2, 2*max(X), 1000)
        plt.plot(xp,p(xp),'-',c='c',lw=1)
        
        
        Y = np.array(SSADListY)
        X = np.array(SSADListX)
            
        SSADslope, SSADintercept, SSADrval, SSADpval, SSADstderr = sc.stats.linregress(X,Y)
        print metric,': r-squared and slope for SSADs:',SSADrval**2, SSADslope
        
        z = np.polyfit(X,Y,1)
        p = np.poly1d(z)
        xp = np.linspace(min(X)/2, 2*max(X), 1000)
        plt.plot(xp,p(xp),'-',c='orange',lw=1)
        
        if index == 0:
            x = 0.5
            y = 8.1
            y2 = 7.4
        elif index == 1:
            x = 0.75
            y = 14.5
            y2 = 13.0
        elif index == 2:
            x = 0.75
            y = -3.0
            y2 = -3.3
        
        plt.text(x, y, 'slope ='+str(round(RADslope,2))+', ' + r'$R^2$' + '=' +str(round(RADrval**2,2)), fontsize=fs-2, color='DarkCyan')         
        plt.text(x, y2, 'slope ='+str(round(SSADslope,2))+', ' + r'$R^2$' + '=' +str(round(SSADrval**2,2)), fontsize=fs-2, color='DarkOrange')         
    
        
        if index == 0:
            plt.scatter([0],[-1], color = 'SkyBlue', alpha = 1, s=10, linewidths=0.9, edgecolor='DarkCyan',
                        label= '1000 randomly drawn Macrobe RADs')
            
            plt.scatter([0],[-1], color = 'PeachPuff',alpha= 1, s=10, linewidths=0.9, edgecolor='DarkOrange',
                        label= '1000 randomly drawn Macrobe SSADs')
            
            plt.legend(bbox_to_anchor=(-0.03, 1.1, 3.9, .2), loc=10, ncol=2,
                                mode="expand",prop={'size':fs})
        
        
        if index > 0: plt.xlabel('Abundance of a sample or an OTU,log', fontsize=fs-2)
        else: plt.xlabel('Number of species or sites, log', fontsize=fs-2)
        
        plt.ylabel(metric, fontsize=fs-2)
        plt.tick_params(axis='both', which='major', labelsize=fs-3)
        
    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    #plt.savefig(mydir+'/figs/Locey_Lennon_Fig2-'+OrC+'_REF.png', dpi=600, bbox_inches = "tight")
    plt.savefig(mydir+'/figs/Locey_Lennon_Fig2-Macro.png', dpi=600, bbox_inches = "tight")
    plt.close()
    
    return




def Fig2Alt():

    Nlist, Slist, Evarlist, ESimplist, ENeelist, EHeiplist, EQlist = [[], [], [],
                                                                    [], [], [], []]
    klist, Shanlist, BPlist, SimpDomlist, SinglesList, tenlist, onelist = [[], [],
                                                                [], [], [], [], []]
    
    NmaxList, rareOnesList, rareRelList, rareSkews = [[], [], [], []]
    
    ct = 0
    klist = []
    
    OrC = 'closed'
    
    SsadDATA = open(mydir+'output/EMP'+OrC+'-SSAD-ResultsTable.txt','r')
    
    for data in SsadDATA:
            
        data_list = data.split()
        N, S, Evar, ESimp, ENee, EHeip, EQ, EPielou, BP, SimpDom, rareRel, rareOnes, skew = data_list
            
        N = float(N)
        S = float(S)
        
        Evarlist.append(float(Evar))    
        Nlist.append(N)
        Slist.append(S)
            
        ESimplist.append(float(ESimp))
        
        BPlist.append(float(BP))
        NmaxList.append(float(BP)*float(N))
        
        EHeiplist.append(float(EHeip))
        rareOnesList.append(float(rareOnes))
        
        rareSkews.append(float(skew))
        
        klist.append('DarkOrange')
        
        ct+=1
        
            
    
    SsadDATA = open(mydir+'output/Macro-SSAD-ResultsTable.txt','r')
    
    
    for data in SsadDATA:
            
        data_list = data.split()
        N, S, Evar, ESimp, ENee, EHeip, EQ, EPielou, BP, SimpDom, rareRel, rareOnes, skew = data_list
            
        N = float(N)
        S = float(S)
        
        Evarlist.append(float(Evar))    
        Nlist.append(N)
        Slist.append(S)
            
        ESimplist.append(float(ESimp))
        
        BPlist.append(float(BP))
        NmaxList.append(float(BP)*float(N))
        
        EHeiplist.append(float(EHeip))
        rareOnesList.append(float(rareOnes))
        
        rareSkews.append(float(skew))
        
        klist.append('DarkCyan')
        
        ct+=1
        

        
    metrics = [['log(number of singletons)', rareOnesList], 
              ['log(greatest abundance)', NmaxList], 
              ['log(evenness)', Evarlist]]
            


    fig = plt.figure()
    for index, i in enumerate(metrics):
    
        ax = fig.add_subplot(2, 3, index+1)
        ax.set_axis_bgcolor('w')
        
        metric = i[0]
        metlist = i[1]
        
        MacSSADListX = []
        MacSSADListY = []
        
        MicSSADListX = []
        MicSSADListY = []
        
        for j, k in enumerate(klist):
            
            if k == 'DarkOrange': 
                if index > 0: MicSSADListX.append(Nlist[j])
                else: MicSSADListX.append(Slist[j])
                
                if metlist[j] == 0: MicSSADListY.append(1)
                else: MicSSADListY.append(metlist[j])
        
                
                
            elif k == 'DarkCyan':
                if index > 0: MacSSADListX.append(Nlist[j])
                else: MacSSADListX.append(Slist[j])
                
                if metlist[j] == 0: MacSSADListY.append(1)
                else: MacSSADListY.append(metlist[j])
                
        
        if index == 0:
            plt.ylim(0, 9)
            plt.xlim(0, 9)
            pass
            
        if index == 1:
            plt.ylim(0, 16)
            plt.xlim(0, 16)
            pass
                
        if index == 2:
            plt.ylim(-3.5, 0.5)
            plt.xlim(0, 16)
            pass
        
        #print len(MacSSADListX), len(MacSSADListY), len(MicSSADListX), len(MicSSADListY)
        #sys.exit()
        
        # scatter plots
        indices = range(1000)
        random.shuffle(indices)            
        
        MacSSADListX = np.log(MacSSADListX)
        MacSSADListY = np.log(MacSSADListY)
        
        MicSSADListX = np.log(MicSSADListX)
        MicSSADListY = np.log(MicSSADListY)
        
        for i in indices:
            
            plt.scatter(MicSSADListX[i], MicSSADListY[i], color = 'PeachPuff', alpha= 1 , s = 4, linewidths=0.5, edgecolor='DarkOrange')                                                
            plt.scatter(MacSSADListX[i], MacSSADListY[i], color = 'SkyBlue', alpha= 1 , s = 4, linewidths=0.5, edgecolor='DarkCyan')
        
        
        Y = np.array(MacSSADListY)
        X = np.array(MacSSADListX)
            
        MacSlope, MacIntercept, MacRval, MacPval, MacStderr = sc.stats.linregress(X,Y)
        print metric,': r-squared and slope for Macrobe SSADs:', MacRval**2, MacSlope
        
        z = np.polyfit(X,Y,1)
        p = np.poly1d(z)
        xp = np.linspace(min(X)/2, 2*max(X), 1000)
        plt.plot(xp,p(xp),'-',c='c',lw=1)
        
        
        Y = np.array(MicSSADListY)
        X = np.array(MicSSADListX)
            
        MicSlope, MicIntercept, MicRval, MicPval, MicStderr = sc.stats.linregress(X,Y)
        print metric,': r-squared and slope for Microbe SSADs:', MicRval**2, MicSlope
        
        z = np.polyfit(X,Y,1)
        p = np.poly1d(z)
        xp = np.linspace(min(X)/2, 2*max(X), 1000)
        plt.plot(xp,p(xp),'-',c='orange',lw=1)
        
        if index == 0:
            x = 0.5
            y = 8.1
            y2 = 7.4
        elif index == 1:
            x = 0.75
            y = 14.5
            y2 = 13.0
        elif index == 2:
            x = 0.75
            y = -3.0
            y2 = -3.3
        
        plt.text(x, y, 'slope ='+str(round(MacSlope,2))+', ' + r'$R^2$' + '=' +str(round(MacRval**2,2)), fontsize=fs-2, color='DarkCyan')         
        plt.text(x, y2, 'slope ='+str(round(MicSlope,2))+', ' + r'$R^2$' + '=' +str(round(MicRval**2,2)), fontsize=fs-2, color='DarkOrange')         
    
        
        if index == 0:
            plt.scatter([0],[-1], color = 'SkyBlue', alpha = 1, s=10, linewidths=0.9, edgecolor='DarkCyan',
                        label= '1000 Macrobe SSADs')
            
            plt.scatter([0],[-1], color = 'PeachPuff',alpha= 1, s=10, linewidths=0.9, edgecolor='DarkOrange',
                        label= '1000 Microbe SSADs')
            
            plt.legend(bbox_to_anchor=(-0.03, 1.1, 3.9, .2), loc=10, ncol=2,
                                mode="expand",prop={'size':fs})
        
        
        if index > 0: plt.xlabel('Abundance of a sample or an OTU, log', fontsize=fs-2)
        else: plt.xlabel('Number of sites, log', fontsize=fs-2)
        
        plt.ylabel(metric, fontsize=fs-2)
        plt.tick_params(axis='both', which='major', labelsize=fs-3)
        
    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    plt.savefig(mydir+'/figs/Locey_Lennon_Fig2Alt.png', dpi=600, bbox_inches = "tight")
    plt.close()
    
    return




def Fig3():

    """ A figure demonstrating a strong abundance relationship across 30
    orders of magnitude in total abundance. The abundance of the most abundant
    species scales in a log-log fashion with the total abundance of the sample
    or system. """

    Nlist, Slist, Evarlist, ESimplist, ENeelist, EHeiplist, EQlist = [[], [], [],
                                                            [], [], [], []]
    klist, Shanlist, BPlist, SimpDomlist, SinglesList, tenlist, onelist = [[], [],
                                                        [], [], [], [], []]
    
    NmaxList, rareOnesList, rareRelList, rarePairList, rareSumOnesList = [[], [],
                                                                    [], [], []]
    
    OrC = 'closed' # is the microbial data (Earth Microbiome Project) going to 
                    # represent closed or open reference OTU assignment

    ct = 0
    klist = []
    
    radDATA = open(mydir+'output/EMP'+OrC+'-RADdata.txt','r')        
    for data in radDATA:
            
        data_list = data.split()
        N, S, ESimp, EHeip, BP, SimpDom, rareRel, rareOnes, rareSumOnes, Var, Evar = data_list
        
        Nlist.append(float(N))
        Slist.append(float(S))
        
        Evarlist.append(float(Evar))   
        NmaxList.append(float(BP)*float(N))    
            
        klist.append('DarkCyan')
        
        ct+=1
        

        
    metrics = [['log'+r'$_{10}$'+'(greatest abundance)', NmaxList], 
              ['log(evenness)', Evarlist]]
            

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
        
        #HGgx = log10((7*10**9)*(10**13))
        #HGgy = log10((7*10**9)*(0.0053*(10**13)))
        
        COWx = log10(2.226*10**15)
        COWy = log10((0.5/80)*(2.226*10**15))
        
        #gCOWx = log10((1.4*10**9) * (2.226*10**15))
        #gCOWy = log10((1.4*10**9) * (0.5/80) * (2.226*10**15))
        
        
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
        
        ax.text(5.5, -3.15, 'log'+r'$_{10}$'+'(Abundance of a sample)', fontsize=fs*1.8)
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
        
    plt.savefig(mydir+'/figs/Locey_Lennon_Fig3v2-'+OrC+'_REF.png', dpi=600,
                bbox_inches = "tight")
    plt.close()
    
    return
    
    
    
    
    

def TL():
    
    """ A figure demonstrating mean-variance relationships for abundances among 
    taxa in the same community (sample) and for abundances of individual taxa
    across samples. Sample data are taken from the Earth Microbiome Project """ 
    
    fig = plt.figure()
    
    Nlist = [] 
    Slist = []
    NmaxList = []
    Varlist = []
    
    ct = 0
    radDATA = open(mydir+'output/EMP'+OrC+'-RADdata.txt','r')
            
    for data in radDATA:
            
        data_list = data.split()
        N, S, ESimp, EHeip, BP, SimpDom, rareRel, rareOnes, rareSumOnes, Var = data_list
        
        N = float(N)
        S = float(S)
        
        Nlist.append(N/S)
        Slist.append(S)
        
        Varlist.append(float(Var))
        NmaxList.append(float(BP)*float(N))
        
        ct+=1
        
            
    ylabel = 'Variance in abundance, log'+r'$_{10}$'
    fig.add_subplot(2, 2, 1)
        
    Varlist = np.log10(Varlist).tolist()
    Nlist = np.log10(Nlist).tolist()
            
    X = []
    Y = []
    for i, val in enumerate(Varlist):
        if np.isnan(val) or np.isinf(val): pass
        else:
            X.append(Nlist[i])
            Y.append(val)    
            
    slope, intercept, rval, pval, stderr = sc.stats.linregress(X,Y)
    print 'r-squared and slope for scatterNmax:',rval**2, slope, intercept
    
    
    plt.hexbin(Nlist, Varlist, mincnt=1, gridsize = 30, bins='log',
             cmap=plt.cm.jet, alpha=1)
    
    plt.text(0.15, 7.8, 'Representing '+str(len(Nlist))+' distributions of\nabundance among OTUs', fontsize=fs-3)         
    
    z = np.polyfit(X,Y,1)
    p = np.poly1d(z)
    xmin =  min(min(X), min(Y))
    xmax =  max(max(X), max(Y))
    
    xp = np.linspace(xmin, xmax, 1000)
    plt.plot(xp, p(xp),'-', c='k', lw=2)
    
    plt.ylim(0, 9)
    plt.xlim(0.0, 2.5)
                
    plt.ylabel(ylabel, fontsize=fs)
    plt.xlabel('Average abundance of OTUs in\nthe same EMP sample, log'+ r'$_{10}$', fontsize=fs)
    
    plt.tick_params(axis='both', which='major', labelsize=fs-1)
                                                            
    plt.text(1.5, 1, 'slope = '+str(round(slope,2))+'\n' + r'$R^2$' + '=' +str(round(rval**2,2)), fontsize=fs)         
      
    ssadDATA = open(mydir+'output/EMP'+OrC+'-SSAD-ResultsTable.txt','r')
                                                          
    Nlist, Slist, Evarlist, ESimplist, ENeelist, EHeiplist, EQlist = [[], [], [],
                                                                    [], [], [], []]
    klist, Shanlist, BPlist, SimpDomlist, SinglesList, tenlist, onelist = [[], [],
                                                                [], [], [], [], []]
    
    NmaxList, rareSkews, rareOnesList = [[], [], []]
    
    Varlist = []
    ct = 0
    
    for data in ssadDATA:
            
        data_list = data.split()
        N, S, Evar, ESimp, ENee, EHeip, EQ, EPielou, BP, SimpDom, rareRel, rareOnes, skew = data_list
            
        N = float(N)
        S = float(S)
        
        Evarlist.append(float(Evar))    
        Nlist.append(N)
        Slist.append(S)
            
        ESimplist.append(float(ESimp))
        
        BPlist.append(float(BP))
        NmaxList.append(float(BP)*float(N))
        
        EHeiplist.append(float(EHeip))
        rareOnesList.append(float(rareOnes))
        
        rareSkews.append(float(skew))
        
        
        ct+=1
        
            
    ylabel = 'Variance in abundance, log'+r'$_{10}$'
    fig.add_subplot(2, 2, 2)
        
    Varlist = np.log10(Varlist).tolist()
    Nlist = np.log10(Nlist).tolist()
        
    X = []
    Y = []
        
    for i, val in enumerate(Varlist):
        if np.isnan(val) or np.isinf(val):
            pass
        else:
            X.append(Nlist[i])
            Y.append(val)    
            
    
    slope, intercept, rval, pval, stderr = sc.stats.linregress(X,Y)
    print 'r-squared and slope for scatterNmax:',rval**2, slope, intercept
    
                 
    plt.ylim(-2, 10)
    plt.xlim(0, 4)
        
    plt.hexbin(Nlist, Varlist, mincnt=1, gridsize = 30, bins='log',
             cmap=plt.cm.jet, alpha=1)
    
    
    plt.text(0.15, 8.3, 'Representing '+str(len(Nlist))+' distributions of\nOTU abundances', fontsize=fs-3)         
    
    
    z = np.polyfit(X,Y,1)
    p = np.poly1d(z)
    xmin =  min(min(X), min(Y))
    xmax =  max(max(X), max(Y))
    
    xp = np.linspace(xmin, xmax, 1000)
    plt.plot(xp, p(xp),'-', c='k', lw=2)
    
    plt.ylabel(ylabel, fontsize=fs)
    plt.xlabel('Average abundance of an OTU\nacross EMP samples, log'+ r'$_{10}$', fontsize=fs)
    
    plt.tick_params(axis='both', which='major', labelsize=fs-1)

    text = 'Left: Samples with greater average OTU abundance have greater variance.\n' 
    text += 'Right: OTUs with greater average sample abundance also have greater variance.\n'
    text += 'Despite having similar slopes, these otherwise different patterns are exceptional\n'
    text += 'to Taylor\'s Law, which predicts a slope between 1 and 2.' 
    
    plt.text(-5.3, 11.5, text, fontsize=fs+1)
    
    plt.text(2.75, -0.8, 'slope = '+str(round(slope,2))+'\n' + r'$R^2$' + '=' +str(round(rval**2,2)), fontsize=fs)
    
    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    plt.savefig(mydir+'/figs/Fig2/SSAD-ScatterMetrics-Heat-'+OrC+'.png', dpi=600, bbox_inches = "tight")
    plt.close()
    
    
    return



def MultipleRegression(KindList, rareSkews, Evarlist, NmaxList, Nlist, Slist, NSlist, headers):
    
    # headers = [ "Category", "Skewnness", "Evenness", "Dominance", "N", "S", "AvgAb"] 
    
    OUT = open(mydir + 'output/MultipleRegressionData.csv','w+')
    
    print >>OUT, "" ","+str(headers[0])+','+str(headers[1])+','+str(headers[2])+','+str(headers[3])+','+str(headers[4])+','+str(headers[5])+','+str(headers[6])
    
    for i, val in enumerate(KindList):
        print >>OUT, str(i+1)+','+str(val)+','+str(rareSkews[i])+','+str(Evarlist[i])+','+str(NmaxList[i])+','+str(Nlist[i])+','+str(Slist[i])+','+str(NSlist[i])
        
    OUT.close()
        
    data = pandas.read_csv(mydir + 'output/MultipleRegressionData.csv', sep=',') #, na_values="."
    
    PanString = str(headers[1])+' ~ '+str(headers[4])
    model = ols(PanString, data).fit()
    print(model.summary())
    
    sys.exit()
    
    return
    

""" The following lines call figure functions to reproduce figures from the 
    Locey and Lennon (2014) manuscript """

Fig1()

#Fig2()
#Fig2Alt()

#Fig3()
#TL()