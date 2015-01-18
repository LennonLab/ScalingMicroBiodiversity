from __future__ import division
#from __future__ import print_function
import  matplotlib.pyplot as plt

import numpy as np
from scipy.stats import gaussian_kde
import random
import scipy as sc
from scipy import stats

import os
import sys
from scipy.stats.distributions import t

import statsmodels.stats.api as sms
import statsmodels.formula.api as smf
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from statsmodels.stats.outliers_influence import summary_table

#import statsmodels.tsa.api as smTsa
import pandas as pd
#import patsy
from math import log10
import linecache



mydir = os.path.expanduser("~/Desktop/Repos/rare-bio/")
mydir2 = os.path.expanduser("~/Desktop/")

sys.path.append(mydir + "tools/feasible_functions")
import feasible_functions as ff

sys.path.append(mydir + "tools/partitions")
#import partitions


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



def Fig1():

    fs = 10 # font size used across figures
    color = str()
    OrC = 'open'
    
    Nlist, Slist, Evarlist, ESimplist, ENeelist, EHeiplist, EQlist = [[], [], [], [], [], [], []]
    klist, Shanlist, BPlist, SimpDomlist, SinglesList, tenlist, onelist = [[], [], [], [], [], [], []]
    NmaxList, rareOnesList, rareRelList, rarePairList, rareSkews, KindList = [[], [], [], [], [], []]
    NSlist = []
    
    ct = 0
    klist = []
    
    radDATA = []
    datasets = []
    BadNames = ['.DS_Store', 'EMPclosed', 'BCI', 'AGSOIL', 'SLUDGE']
            
    for name in os.listdir(mydir2 +'data/micro'):
            if name in BadNames: continue
            
            path = mydir2+'data/micro/'+name+'/'+name+'-SADMetricData.txt'
            num_lines = sum(1 for line in open(path))
            datasets.append([name, 'micro', num_lines])   
    
    for name in os.listdir(mydir2 +'data/macro'):
            if name in BadNames: continue
            
            path = mydir2+'data/macro/'+name+'/'+name+'-SADMetricData.txt'
            num_lines = sum(1 for line in open(path))
            datasets.append([name, 'macro', num_lines])   
    
    numMac = 0
    numMic = 0
            
    for dataset in datasets:
        
        name, kind, numlines = dataset
        lines = []
        
        if numlines > 40: lines = random.sample(range(1, numlines+1), 40)
        else: lines = random.sample(range(1, numlines+1), 40)
        
        path = mydir2+'data/'+kind+'/'+name+'/'+name+'-SADMetricData.txt'
        
        radDATA = []
        for line in lines:
            data = linecache.getline(path, line)
            radDATA.append(data)
        
        print name, kind, numlines, len(radDATA)
        
        for data in radDATA:
            
            data = data.split()
            if len(data) == 0: 
                print 'no data'
                continue
            
            name, kind, N, S, Evar, ESimp, ENee, EHeip, EQ, EPielou, BP, SimpDom, rareRel, rareOnes, skew = data
            
            N = float(N)
            S = float(S)
            
            if S < 2: continue # Min species richness
            
            Nlist.append(float(np.log(N)))
            Slist.append(float(np.log(S)))
            NSlist.append(float(np.log(N/S)))
            
            Evarlist.append(float(np.log(float(Evar))))
            ESimplist.append(float(np.log(float(ESimp))))
            KindList.append(kind)
            
            BPlist.append(float(BP))
            NmaxList.append(float(np.log(float(BP)*float(N))))
            EHeiplist.append(float(EHeip))
            
            rareOnesList.append(float(rareOnes))
            rareRelList.append(float(rareOnes)/S)
            
            # lines for the log-modulo transformation of skewnness
            skew = float(skew)
            sign = 1
            if skew < 0: sign = -1
            
            lms = np.log(np.abs(skew) + 1)
            lms = lms * sign
            #if lms > 3: print name, N, S
            rareSkews.append(float(lms))
            
            if kind == 'micro':
                numMic += 1 
                klist.append('b')
            if kind == 'macro': 
                klist.append('r')
                numMac += 1
        
            ct+=1
    
    print 'Mic:',numMic,'Mac:', numMac
    
    
    metrics = [['Rarity, log-modulo(skewnness)', rareSkews], 
              ['Dominance, log(greatest abundance)', NmaxList], 
              ['Evenness, log(Simpsons)', ESimplist]]
            

    fig = plt.figure()
    for index, i in enumerate(metrics):
    
        fig.add_subplot(2, 3, index+1)
        
        metric = i[0]
        metlist = i[1]
        
        # EXAMINE OLS RESULTS BY SHUFFLING THESE THREE LISTS     
        #np.random.shuffle(KindList)
        #np.random.shuffle(Nlist)
        #np.random.shuffle(metlist)
        
        MacListX = []
        MacListY = []
        
        MicListX = []
        MicListY = []
        
        nmacs = 0
        nmics = 0
        
        for j, k in enumerate(KindList):
            
            if k == 'micro': 
                nmics += 1
                if index >= 0: MicListX.append(Nlist[j])
                else: MicListX.append(Slist[j])
                
                MicListY.append(metlist[j])
        
            elif k == 'macro':
                nmacs += 1
                if index >= 0: MacListX.append(Nlist[j])
                else: MacListX.append(Slist[j])
                
                MacListY.append(metlist[j])
        
        if index == 0:
            plt.ylim(0, 5)
            plt.xlim(0, 16)
            pass
            
        if index == 1:
            plt.ylim(0, 14)
            plt.xlim(0, 16)
            pass
                
        if index == 2:
            plt.ylim(-3.5, 0.1)
            plt.xlim(0, 16)
            pass
        
        micY = []
        micX = []
        macY = []
        macX = []
    
        for i in range(len(MacListX)):
            plt.scatter(MacListX[i], MacListY[i], color = 'LightCoral', alpha= 1 , s = 4, linewidths=0.5, edgecolor='Crimson')                                                
            plt.scatter(MicListX[i], MicListY[i], color = 'SkyBlue', alpha= 1 , s = 4, linewidths=0.5, edgecolor='Steelblue')
            
            micY.append(MicListY[i])    
            macY.append(MacListY[i])    
            micX.append(MicListX[i])    
            macX.append(MacListX[i]) 
        
        ### OLS for macrobes
        #Y = np.array(macY+micY)
        #X = np.array(macX+micX)
        #MacSlope, MacInt, MacRval, MacPval, MacStderr = sc.stats.linregress(X,Y)
        #print metric,' for Combined: r-squared, slope, and p-val for macrobes:',MacRval**2, MacSlope, MacPval
        
        # Multiple regression
        d = pd.DataFrame({'N': list(Nlist)})
        d['y'] = list(metlist)
        d['Kind'] = list(KindList)

        f = smf.ols('y ~ N * Kind', d).fit()
        
        prstd, iv_l, iv_u = wls_prediction_std(f)
        
        MacPIx = []
        MacPIl = []
        MacPIh = []
        MacFitted = []
        
        MicPIx = []
        MicPIl = []
        MicPIh = []
        MicFitted = []
        
        for j, kval in enumerate(KindList):
            if kval == 'macro':
                MacPIx.append(Nlist[j])
                MacPIl.append(iv_l[j])
                MacPIh.append(iv_u[j])
                MacFitted.append(f.fittedvalues[j])
                
            elif kval == 'micro':
                MicPIx.append(Nlist[j])
                MicPIl.append(iv_l[j])
                MicPIh.append(iv_u[j])
                MicFitted.append(f.fittedvalues[j])
        
        plt.plot(MacPIx, MacPIh, ls=':', lw=0.1, c='r', alpha=0.6)
        plt.plot(MacPIx, MacPIl, ls=':', lw=0.1, c='r', alpha=0.6)
        #plt.fill_between(MacPIx, MacPIl, MacPIh, linewidth=0.0, facecolor= 'r', alpha=0.9)
        plt.plot(MacPIx, MacFitted, ls='-', c='r', alpha=0.6)
        
        plt.plot(MicPIx, MicPIh, ls=':', lw=0.1, c='b', alpha=0.6)
        plt.plot(MicPIx, MicPIl, ls=':', lw=0.1, c='b', alpha=0.6)
        #plt.fill_between(MicPIx, MicPIl, MicPIh, linewidth=0.0, facecolor= 'b', alpha=0.9)
        plt.plot(MicPIx, MicFitted, ls='-', c='b', alpha=0.6)
        
        print '\n',metric
        print f.summary(), '\n'
        
        """
        HC = sms.linear_harvey_collier(f) # Harvey Collier test for linearity. The Null hypothesis is that the regression is correctly modeled as linear.
        print 'Harvey-Collier test for linearity:', HC
        
        RB = sms.linear_rainbow(f) # Rainbow test for linearity. The Null hypothesis is that the regression is correctly modeled as linear.
        print 'Rainbow test for linearity:', RB
        
        LM = sms.linear_lm(f.resid, f.model.exog)
        print 'Lagrangian multiplier test for linearity:', LM
        
        BGtest = sms.acorr_breush_godfrey(f, nlags=None, store=False) # Breusch Godfrey Lagrange Multiplier tests for residual autocorrelation
        print 'Breusch Godfrey test for autocorrelation:', BGtest # Lagrange multiplier test statistic, p-value for Lagrange multiplier test, fstatistic for F test, pvalue for F test
        
        NeweyWest = sms.sandwich_covariance.cov_hac(f)
        print NeweyWest
        """
        
        if index == 0:
            x = 2
            x2 = 2
            y = 5.6
            y2 = 5.2
            y3 = 4.1
            
        elif index == 1:
            x = 2
            x2 = 2
            y = 15.7
            y2 = 14.6
            y3 = 11.8
            
        elif index == 2:
            x = 2
            x2 = 9
            y = 0.51
            y2 = 0.26
            y3 = -0.5
        
        params = f.params
        Const = params[0]
        MicroCoef = params[1]
        N_Coef = params[2]
        IntCoef = params[3] 
        
        r2 = f.rsquared
        
        pvals = f.pvalues
        ConstPval = pvals[0]
        MicroPval = pvals[1]
        N_Pval = pvals[2]
        IntPval = pvals[3]
        
        if index == 0:
            plt.text(x, y, r'$y_{micro}$'+ ' = '+str(round(Const+MicroCoef,2))+'+'+str(round(N_Coef, 2))+'*'+r'$N$', fontsize=fs-2, color='Steelblue')         
            plt.text(x, y2, r'$y_{macro}$'+ ' = '+str(round(Const,2))+'+'+str(round(N_Coef,2))+'*'+r'$N$', fontsize=fs-2, color='Crimson')  
        
        elif index == 1:
            plt.text(x, y, r'$y_{micro}$'+ ' = '+str(round(Const,2))+'+'+str(round(N_Coef,2))+'*'+r'$N$', fontsize=fs-2, color='Steelblue')         
            plt.text(x, y2, r'$y_{macro}$'+ ' = '+str(round(Const,2))+'+'+str(round(N_Coef,2))+'*'+r'$N$', fontsize=fs-2, color='Crimson')         
        
        elif index == 2:
            plt.text(x, y, r'$y_{micro}$'+ ' = '+str(round(Const+MicroCoef,2))+' '+str(round(N_Coef+IntCoef,2))+'*'+r'$N$', fontsize=fs-2, color='Steelblue')         
            plt.text(x,y2, r'$y_{macro}$'+ ' = '+str(round(Const,2))+' '+str(round(N_Coef,2))+'*'+r'$N$', fontsize=fs-2, color='Crimson')         
        
        plt.text(x2, y3,  r'$R^2$' + '=' +str(round(r2,3)), fontsize=fs-2, color='k')  
        
        if index == 0:
            plt.scatter([0],[-1], color = 'SkyBlue', alpha = 1, s=10, linewidths=0.9, edgecolor='Steelblue',
                        label= 'microbes (n='+str(len(MicListY))+')')
            
            plt.scatter([0],[-1], color = 'LightCoral',alpha= 1, s=10, linewidths=0.9, edgecolor='Crimson',
                        label= 'macrobes (n='+str(len(MacListY))+')')
            
            plt.legend(bbox_to_anchor=(-0.04, 1.2, 3.89, .2), loc=10, ncol=2,
                                mode="expand",prop={'size':fs})
        
        
        plt.xlabel('log(total abundance of a sample)', fontsize=fs-2)
        
        plt.ylabel(metric, fontsize=fs-2)
        plt.tick_params(axis='both', which='major', labelsize=fs-3)
        
    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    plt.savefig(mydir+'/figs/Locey_Lennon_Fig1-'+OrC+'_REF.png', dpi=600, bbox_inches = "tight")
    plt.close()
    
    return






def Fig2():

    fs = 10 # font size used across figures
    color = str()
    OrC = 'open'

    Nlist, Slist, Evarlist, ESimplist, ENeelist, EHeiplist, EQlist = [[], [], [],
                                                                    [], [], [], []]
    klist, Shanlist, BPlist, SimpDomlist, SinglesList, tenlist, onelist = [[], [],
                                                                [], [], [], [], []]
    
    NmaxList, rareOnesList, rareRelList, rareSkews, NSlist = [[], [], [], [], []]
    
    ct = 0
    klist = []
    
    radDATA = []
    datasets = []
    BadNames = ['.DS_Store', 'EMPclosed', 'BCI', 'AGSOIL', 'SLUDGE', 'CHU', 'HYDRO', 'CATLIN', 'FECES', 'FUNGI', 'LAUB']
            
    for name in os.listdir(mydir2 +'data/micro'):
            if name in BadNames: continue
            
            path = mydir2+'data/micro/'+name+'/'+name+'-SSADMetricData.txt'
            num_lines = 5594412
            #num_lines = sum(1 for line in open(path))
            datasets.append([name, 'micro', num_lines])   
    
    for name in os.listdir(mydir2 +'data/macro'):
            if name in BadNames: continue
            
            path = mydir2+'data/macro/'+name+'/'+name+'-SSADMetricData.txt'
            num_lines = sum(1 for line in open(path))
            datasets.append([name, 'macro', num_lines])   
    
    numMac = 0
    numMic = 0
    KindList = []    
    
    for dataset in datasets:
        
        name = dataset[0]
        kind = dataset[1]
        numlines = dataset[2]    
        
        lines = []
        
        SampSize = 1800
        
        if kind == 'macro' and numlines > 200: lines = random.sample(range(1, numlines+1), 300)
        elif kind == 'macro': lines = random.sample(range(1, numlines+1), 40)
        elif kind == 'micro': lines = random.sample(range(1, numlines+1), SampSize)
        path = mydir2+'data/'+kind+'/'+name+'/'+name+'-SSADMetricData.txt'
        
        radDATA = []
        for line in lines:
            data = linecache.getline(path, line)
            radDATA.append(data)
        
        print name, kind, numlines, len(radDATA)
        
        for data in radDATA:
            
            data = data.split()
            if len(data) == 0: 
                print 'no data'
                continue
            
            N, S, Evar, ESimp, ENee, EHeip, EQ, EPielou, BP, SimpDom, rareRel, rareOnes, skew = data
            
            N = float(N)
            S = float(S)
        
            if S < 1: continue
        
            Evarlist.append(np.log(float(Evar)))    
            Nlist.append(float(np.log(N)))
            Slist.append(float(np.log(S)))
            NSlist.append(float(np.log(N/S)))
            
            ESimplist.append(float(ESimp))
        
            BPlist.append(float(BP))
            NmaxList.append(float(np.log(float(BP)*float(N))))
        
            EHeiplist.append(float(EHeip))
            rareOnesList.append(float(rareOnes))
        
            # log-modulo of skew
            skew = float(skew)
            rareSkews.append(skew)
            
            #sign = 1
            #if skew < 0: sign = -1
            
            #lms = np.log(np.abs(skew) + 1)
            #lms = lms * sign
            #rareSkews.append(float(lms))
        
            if kind == 'micro': 
                klist.append('DarkOrange')
                numMic += 1
            elif kind == 'macro': 
                klist.append('DarkCyan')
                numMac += 1
            
            KindList.append(kind)
            ct+=1
        
    
    print 'Macrobes', numMac, 'Microbes:', numMic
         
    metrics = [['Rarity, log-modulo(skewnness)', 'Rarity', rareSkews], 
              ['Dominance, log(greatest abundance)', 'Dominance', NmaxList], 
              ['Evenness, Simpsons', 'Evenness', ESimplist]]    
    
    fig = plt.figure()
    for i, val in enumerate(metrics):
    
        ax = fig.add_subplot(2, 3, i+1)
        ax.set_axis_bgcolor('w')
        
        ylabel = val[0]
        metric = str(val[1])
        metlist = val[2]
        
        MacSSADListX = []
        MacSSADListY = []
        
        MicSSADListX = []
        MicSSADListY = []
        
        for j, k in enumerate(klist):
            
            if k == 'DarkOrange': 
                MicSSADListX.append(Nlist[j])
                MicSSADListY.append(metlist[j])
        
            elif k == 'DarkCyan':
                MacSSADListX.append(Nlist[j])
                MacSSADListY.append(metlist[j])
                
        """
        if i == 0:
            plt.ylim(-1, 4)
            plt.xlim(0, 16)
            pass
            
        if i == 1:
            plt.ylim(0, 12)
            plt.xlim(0, 16)
            pass
                
        if i == 2:
            plt.ylim(-3.5, 0.5)
            plt.xlim(0, 16)
            pass
        """
        for ii in range(len(MicSSADListX)):
            plt.scatter(MicSSADListX[ii], MicSSADListY[ii], color = 'PeachPuff', alpha= 1 , s = 4, linewidths=0.25, edgecolor='DarkOrange')                                                
            plt.scatter(MacSSADListX[ii], MacSSADListY[ii], color = 'SkyBlue', alpha= 1 , s = 4, linewidths=0.25, edgecolor='DarkCyan')

        
        d = pd.DataFrame({'N': list(Nlist)})
        d['y'] = list(metlist)
        d['Kind'] = list(KindList)

        f = smf.ols('y ~ N * Kind', d).fit()
        print (f.summary()), '\n\n'
        
        prstd, iv_l, iv_u = wls_prediction_std(f)
        
        MacPIx = []
        MacPIl = []
        MacPIh = []
        MacFitted = []
        
        MicPIx = []
        MicPIl = []
        MicPIh = []
        MicFitted = []
        
        for j, kval in enumerate(KindList):
            if kval == 'macro':
                MacPIx.append(Nlist[j])
                MacPIl.append(iv_l[j])
                MacPIh.append(iv_u[j])
                MacFitted.append(f.fittedvalues[j])
                
            elif kval == 'micro':
                MicPIx.append(Nlist[j])
                MicPIl.append(iv_l[j])
                MicPIh.append(iv_u[j])
                MicFitted.append(f.fittedvalues[j])
        
        
        plt.plot(MacPIx, MacPIh, ls='--', lw=0.1, c='c', alpha=0.9)
        plt.plot(MacPIx, MacPIl, ls='--', lw=0.1, c='c', alpha=0.9)
        #plt.fill_between(MacPIx2, MacPIh2, MacPIl2, facecolor= 'c', alpha=0.3)
        plt.plot(MacPIx, MacFitted, ls='-', lw=1, c='c', alpha=0.9)
        
        plt.plot(MicPIx, MicPIh, ls='--', lw=0.1, c='orange', alpha=0.9)
        plt.plot(MicPIx, MicPIl, ls='--', lw=0.1, c='orange', alpha=0.9)
        #plt.fill_between(MicPIx2, MicPIh2, MicPIl2, facecolor= 'orange', alpha=0.3)
        plt.plot(MicPIx, MicFitted, ls='-', lw=1, c='orange', alpha=0.9)
        
        """
        if i == 0:
            x = 1
            y2 = 4.2
            y3 = 3.32
            
        elif i == 1:
            x = 1
            y2 = 12.4
            y3 = 10.3
            
        elif i == 2:
            x = 9
            y2 = 1.14
            y3 = 0.94
        """
        
        params = f.params
        Const = params[0]
        MicroCoef = params[1]
        N_Coef = params[2]
        IntCoef = params[3] 
        
        r2 = f.rsquared
        
        pvals = f.pvalues
        ConstPval = pvals[0]
        MicroPval = pvals[1]
        N_Pval = pvals[2]
        IntPval = pvals[3]
        
        """
        if metric == 'Dominance':
            plt.text(x, y2, r'$y_i$'+ ' = '+str(round(Const,2))+' '+str(round(N_Coef,2))+'*'+r'$N$', fontsize=fs-2, color='DarkCyan')         
            plt.text(x, y2, r'$y_i$'+ ' = '+str(round(Const,2))+' '+str(round(N_Coef,2))+'*'+r'$N$', fontsize=fs-2, color='DarkOrange')
        elif metric == 'Rarity':
            plt.text(x, y2, r'$y_i$'+ ' = '+str(round(Const,2))+' '+str(round(N_Coef,2))+'*'+r'$N$', fontsize=fs-2, color='DarkCyan')         
            plt.text(x, y2, r'$y_i$'+ ' = '+str(round(Const,2))+' '+str(round(N_Coef,2))+'*'+r'$N$', fontsize=fs-2, color='DarkOrange')
        elif metric == 'Evenness':
            plt.text(x, y2, r'$y_i$'+ ' = '+str(round(Const,2))+' '+str(round(N_Coef,2))+'*'+r'$N$', fontsize=fs-2, color='DarkCyan')         
            plt.text(x, y2, r'$y_i$'+ ' = '+str(round(Const,2))+' '+str(round(N_Coef,2))+'*'+r'$N$', fontsize=fs-2, color='DarkOrange')
        """
        #plt.text(x, y2, r'$y_i$'+ '='+str(Const)+' + '+str(N_Coef)+'*'+r'$N$'+' + '+str(MicroCoef)+'*'+'r$micro$', fontsize=fs-3, color='k')         
        
        #plt.text(x, y3,  r'$R^2$' + '=' +str(round(r2,3)), fontsize=fs-2, color='k')         
        
        """
        if i == 0:
            plt.scatter([0],[-1], color = 'SkyBlue', alpha = 1, s=10, linewidths=0.9, edgecolor='DarkCyan',
                        label= str(SampSize)+' Macrobe SSADs')
            
            plt.scatter([0],[-1], color = 'PeachPuff',alpha= 1, s=10, linewidths=0.9, edgecolor='DarkOrange',
                        label= str(SampSize)+' Microbe SSADs')
            
            plt.legend(bbox_to_anchor=(-0.03, 1.2, 3.9, .2), loc=10, ncol=2,
                                mode="expand",prop={'size':fs})
        """
        plt.xlabel('log(abundance of a taxa in a dataset)', fontsize=fs-3)
        
        plt.ylabel(ylabel, fontsize=fs-3)
        plt.tick_params(axis='both', which='major', labelsize=fs-3)
        
    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    plt.savefig(mydir+'/figs/Locey_Lennon_'+str(OrC)+'_Fig2.png', dpi=600, bbox_inches = "tight")
    #plt.show()
    
    return




def Fig3():

    """ A figure demonstrating a strong abundance relationship across 30
    orders of magnitude in total abundance. The abundance of the most abundant
    species scales in a log-log fashion with the total abundance of the sample
    or system. """

    fs = 10 # font size used across figures
    color = str()

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
    
    fs = 10 # font size used across figures
    color = str()
    
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



""" The following lines call figure functions to reproduce figures from the 
    Locey and Lennon (2014) manuscript """

#Fig1()
Fig2()
#Fig3()
#TL()