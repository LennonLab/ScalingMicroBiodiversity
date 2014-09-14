from __future__ import division
import  matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
import sys
import os
import feasible_functions as ff
import scipy as sc
from scipy.stats import gaussian_kde, sem
import random

from matplotlib import rcParams


def get_kdens(summands):
    """ Finds the kernel density function across a sample of parts
    """
    density = gaussian_kde(summands)
    n = len(summands)
    xs = np.linspace(float(min(summands)),float(max(summands)),n)
    density.covariance_factor = lambda : .5
    density._compute_covariance()
    D = [xs, density(xs)]
    return D



radDATA = open('/Users/lisalocey/Desktop/RareBio/RADdataEMP.txt','r')

fs = 30
color = str()

Nlist = []
Slist = []
Evarlist = [] 
ESimplist = []
ENeelist = []
EHeiplist = []
EQlist = []
klist = []
Shanlist = []

BPlist = []
SimpDomlist = []

SinglesList = []
atylist = []
tenlist = []
onelist = []
ptonelist = []
ptzonelist = []

BP1list = []
BP2list = []
BP3list = []
BP4list = []
BP5list = []
NmaxList = []

Nper90 = []
Nper70 = []
Nper50 = []
Nper30 = []
Nper10 = []


ct = 0
klist = []
for data in radDATA:
    
    data_list = data.split()
    name, kind, N, S, Evar, ESimp, ENee, EHeip, EQ, EPielou, BP, SimpDom, perc_ones, s_80, s_ten, s_one, s_ptone, s_ptzone, BP1, BP2, BP3, BP4, BP5, Nmax, N90, N70, N50, N30, N10 = data_list
    
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
    atylist.append(float(s_80))
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
    
    Nper90.append(float(N90))
    Nper70.append(float(N70))
    Nper50.append(float(N50))
    Nper30.append(float(N30))
    Nper10.append(float(N10))
    
    if kind == 'micro':
        color = 'Yellow'
    elif kind == 'macro':
        color = 'Cyan'
    
    klist.append(color)
    ct+=1
    
Npercs = [NmaxList, Nper90, Nper70, Nper50, Nper30, Nper10]



def kdensMetrics():
    
    metrics = [['Evar', Evarlist], ['Evar', Evarlist], ['Simpsons', ESimplist], ['Heips', EHeiplist], ['Pielous', Shanlist]]
    metrics.extend([['Berger-Parker', BPlist],['Simpson\'s D', SimpDomlist]])

    for index, i in enumerate(metrics):
        
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        
        ax.spines['top'].set_color('w')
        ax.spines['bottom'].set_color('w')
        ax.spines['left'].set_color('w')
        ax.spines['right'].set_color('w')        
        ax.xaxis.label.set_color('w')
        ax.tick_params(axis='both', colors='w')            
        
        metric = i[0]
        metlist = i[1]

        microlist = []
        macrolist = []
        
        for i, j in enumerate(metlist):
            if klist[i] == 'Yellow':
                microlist.append(j)
            elif klist[i] == 'Cyan':
                macrolist.append(j)
        
        D = get_kdens(microlist)
        n = len(microlist)
        plt.plot(D[0], D[1], color = 'Yellow', lw=3)
        ymax1 = max(D[1])
        
        plt.scatter([0],[-1], color = 'Yellow', s=10,
                     label= 'microbial (n='+str(n)+')', alpha=0.0)
        
        
        n = len(macrolist)
        D = get_kdens(macrolist)
        plt.plot(D[0], D[1], color = 'Cyan', lw=3)
        ymax2 = max(D[1])
        
        plt.scatter([0],[-1], color = 'Cyan', s=10,
                     label= 'macrobial (n='+str(n)+')', alpha=0.0)
        
        l = plt.legend(loc=2, ncol=1,
                   prop={'size':fs-5})    
                
        l.get_frame().set_alpha(0.0)
        txtclrs = ['Yellow', 'Cyan']
        for ti, text in enumerate(l.get_texts()):
            text.set_color(txtclrs[ti])
        
        
        ymax = max(ymax1, ymax2) + 0.4*(max(ymax1, ymax2))
        plt.ylim(0.0, ymax)
        plt.ylabel("Probability density", fontsize=fs, color='w')
        plt.xlabel(metric, fontsize=fs, color='w')
        plt.tick_params(axis='both', which='major', labelsize=fs-10)

        plt.subplots_adjust(wspace=0.4, hspace=0.4)
        
        ax.spines['top'].set_color('w')
        ax.spines['bottom'].set_color('w')
        ax.spines['left'].set_color('w')
        ax.spines['right'].set_color('w')        
        ax.xaxis.label.set_color('w')
        ax.tick_params(axis='both', colors='w', length=8, width=4)
        
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(4)
        
        rcParams.update({'figure.autolayout': True})
        plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/'+metric+'kdens.png',
        transparent=True, dpi=600, pad_inches = 0.1, layout='tight')
        plt.close()
        print 'kdens fig for '+metric+': done'   
    
    return
    
    
def histMetrics():
    
    metrics = [['Evar', Evarlist], ['Evar', Evarlist], ['Simpsons', ESimplist], ['Heips', EHeiplist], ['Pielous', Shanlist]]
    metrics.extend([['Berger-Parker', BPlist],['Simpson\'s D', SimpDomlist]])

    for index, i in enumerate(metrics):
        fig = plt.figure()
    
        ax = fig.add_subplot(1, 1, 1)
        
        ax.spines['top'].set_color('w')
        ax.spines['bottom'].set_color('w')
        ax.spines['left'].set_color('w')
        ax.spines['right'].set_color('w')        
        ax.xaxis.label.set_color('w')
        ax.tick_params(axis='both', colors='w', length=8, width=4)
        
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(4)
        
        metric = i[0]
        metlist = i[1]

        microlist = []
        macrolist = []
        
        for i, j in enumerate(metlist):
            if klist[i] == 'Yellow':
                microlist.append(j)
            elif klist[i] == 'Cyan':
                macrolist.append(j)
        
        plt.hist(macrolist, 40, color='Cyan', linewidth='3', histtype='step',
        label = 'macrobial (n='+str(len(macrolist))+')')
        
        plt.hist(microlist, 40, color='Yellow', linewidth='3', histtype='step',
        label = 'microbial (n='+str(len(microlist))+')')
        
        #l = plt.legend(bbox_to_anchor=(-0.015, 1, 1.029, .2), loc=1, ncol=2,
        #                       mode="expand", prop={'size':fs-5}) 
        
        l = plt.legend(loc=1, ncol=1,
                   prop={'size':fs-5})     
                    
        l.get_frame().set_alpha(0.0)
        for text in l.get_texts():
            text.set_color('w')
        
        
        #ymax = max(ymax1, ymax2) + 0.4*(max(ymax1, ymax2))
        #plt.ylim(0.0, ymax)
        plt.xlim(0.0, 1.0)
        plt.ylabel("Number of communities",fontsize=fs, color='w')
        plt.xlabel(metric+' evenness', fontsize=fs, color='w')
        plt.tick_params(axis='both', which='major', labelsize=15)

        #rcParams.update({'figure.autolayout': True})
        
        plt.subplots_adjust(wspace=0.4, hspace=0.4)
        
        ax.spines['top'].set_color('w')
        ax.spines['bottom'].set_color('w')
        ax.spines['left'].set_color('w')
        ax.spines['right'].set_color('w')        
        ax.xaxis.label.set_color('w')
        ax.tick_params(axis='both', colors='w', length=8, width=4)
        
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(4)
        
        rcParams.update({'figure.autolayout': True})
        plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/hist'+metric+'.png',
        transparent=True, dpi=600, pad_inches = 0.3)
        
        plt.close()
        print 'hist fig for '+metric+': done'   
    
    return
    
    
    
def scatterMetrics():

    metrics = [['Evar', Evarlist], ['Simpsons', ESimplist], ['Heips', EHeiplist], ['Pielous', Shanlist]]
    metrics.extend([['Berger-Parker', BPlist],['Simpson\'s D', SimpDomlist]])

    Nlocs = [1, 1, 1, 3, 2, 3]
    Slocs = [1, 1, 1, 4, 1, 4]
    
    rcParams.update({'figure.autolayout': True})
        
    for j, i in enumerate(metrics):
        
        nloc = Nlocs[j]
        sloc = Slocs[j]
        
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        
        ax.spines['top'].set_color('w')
        ax.spines['bottom'].set_color('w')
        ax.spines['left'].set_color('w')
        ax.spines['right'].set_color('w')        
        ax.xaxis.label.set_color('w')
        ax.tick_params(axis='both', colors='w', length=8, width=4)
        
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(4)    
        
        metric = i[0]
        metlist = i[1]
        
        nmacs = 0
        nmics = 0
        for j, k in enumerate(klist):
            plt.scatter(Nlist[j], metlist[j], color = k, s = 20,
            linewidths = 0.0)
            if k == 'Yellow': nmics += 1
            elif k == 'Cyan': nmacs += 1
            
        plt.scatter([0],[-1], color = 'Yellow', s=10,
                     label= 'microbial (n='+str(nmics)+')', alpha=0.0)
        plt.scatter([0],[-1], color = 'Cyan', s=10,
                     label= 'macrobial (n='+str(nmacs)+')', alpha=0.0)
        
        l = plt.legend(loc=nloc, ncol=1,
                   prop={'size':fs-10})    
                
        l.get_frame().set_alpha(0.0)
        txtclrs = ['Yellow', 'Cyan']
        for ti, text in enumerate(l.get_texts()):
            text.set_color(txtclrs[ti])
        
        
        plt.ylim(0.0, 1.0)
        plt.xlabel('Total abundance', fontsize=fs, color='w')
        plt.ylabel(metric, fontsize=fs, color='w')
        plt.tick_params(axis='both', which='major', labelsize=fs-10)
        plt.xscale('log')
        
        plt.subplots_adjust(wspace=0.4, hspace=0.4)
        
        plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/'+metric+'_N_scatter.png',
        transparent=True, dpi=600, pad_inches = 0.1)
        plt.close()
          
            
                
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        
        ax.spines['top'].set_color('w')
        ax.spines['bottom'].set_color('w')
        ax.spines['left'].set_color('w')
        ax.spines['right'].set_color('w')        
        ax.xaxis.label.set_color('w')
        ax.tick_params(axis='both', colors='w', length=8, width=4)
        
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(4)
            
        
        for j, k in enumerate(klist):
            plt.scatter(Slist[j], metlist[j], color = k, s = 20,
            linewidths = 0.0)
        
        plt.scatter([0],[-1], color = 'Yellow', s=10,
                     label= 'microbial (n='+str(nmics)+')', alpha=0.0)
        plt.scatter([0],[-1], color = 'Cyan', s=10,
                     label= 'macrobial (n='+str(nmacs)+')', alpha=0.0)
        
        l = plt.legend(loc=nloc, ncol=1,
                   prop={'size':fs-10})    
                
        l.get_frame().set_alpha(0.0)
        txtclrs = ['Yellow', 'Cyan']
        for ti, text in enumerate(l.get_texts()):
            text.set_color(txtclrs[ti])
           
        
        plt.xlabel('Richness', fontsize=fs, color='w')
        plt.ylabel(metric, fontsize=fs, color='w')
        plt.tick_params(axis='both', which='major', labelsize=fs-10)
        plt.xlim(10,5000)
        plt.ylim(0.0, 1.0)
        plt.xscale('log')
    
        plt.subplots_adjust(wspace=0.4, hspace=0.4)
        
        plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/'+metric+'_S_scatter.png',
        transparent=True, dpi=600, pad_inches = 0.1)
        plt.close()
        
        print metric,'scatter plot: done'
    return



def scatterDominance():
    
    metrics = [['Berger-Parker', BPlist], ['Simpsons Dominance', SimpDomlist]]
    #metrics = [['Nmax', NmaxList]]
    
    for i in metrics:
    
        fig = plt.figure()
        ax = fig.add_subplot(2, 2, 1)
        
        ax.spines['top'].set_color('w')
        ax.spines['bottom'].set_color('w')
        ax.spines['left'].set_color('w')
        ax.spines['right'].set_color('w')        
        ax.xaxis.label.set_color('w')
        ax.tick_params(axis='both', colors='w', length=8, width=4)
        
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(4)
        
        metric = i[0]
        metlist = i[1]
        
        nmacs = 0
        nmics = 0
        for j, k in enumerate(klist):
            plt.scatter(Nlist[j], metlist[j], color = k, s = 5,
            linewidths = 0.0)
            if k == 'Yellow': nmics += 1
            elif k == 'Cyan': nmacs += 1
            
        plt.scatter([0],[-1], color = 'Yellow', s=10,
                     label= 'microbial (n='+str(nmics)+')')
        plt.scatter([0],[-1], color = 'Cyan', s=10,
                     label= 'macrobial (n='+str(nmacs)+')')
        
        #l = plt.legend(bbox_to_anchor=(-0.05, 1, 2.5, .2), loc=1, ncol=2,
        #                        mode="expand", prop={'size':fs-5}) 
        
        l = plt.legend(loc=1, ncol=1,
                   prop={'size':fs-5})         
                                
        l.get_frame().set_alpha(0.0)
        for text in l.get_texts():
            text.set_color('w')
        

        plt.ylim(0.0, 1.0)
        plt.xlabel('Total abundance', fontsize=fs, color='w')
        plt.ylabel(metric, fontsize=fs, color='w')
        plt.tick_params(axis='both', which='major', labelsize=fs-3)
        plt.xscale('log')
        
        ax = fig.add_subplot(2, 2, 2)
        
        ax.spines['top'].set_color('w')
        ax.spines['bottom'].set_color('w')
        ax.spines['left'].set_color('w')
        ax.spines['right'].set_color('w')        
        ax.xaxis.label.set_color('w')
        ax.tick_params(axis='both', colors='w', length=8, width=4)
        
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(4)
        
        for j, k in enumerate(klist):
            plt.scatter(Slist[j], metlist[j], color = k, s = 5,
            linewidths = 0.0)


        plt.xlabel('Richness', fontsize=fs, color='w')
        plt.ylabel(metric, fontsize=fs, color='w')
        plt.tick_params(axis='both', which='major', labelsize=fs-3)
        plt.xlim(10,5000)
        plt.ylim(0.0, 1.0)
        plt.xscale('log')
    
        plt.subplots_adjust(wspace=0.4, hspace=0.4)
        
        rcParams.update({'figure.autolayout': True})
        plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/'+metric+'_scatter.png',
        transparent=True, dpi=600, pad_inches = 0.1)
        plt.close()
        
        print metric,'scatter plot: done'
    return


def scatterNmax():
    
    metrics = [['Nmax', NmaxList]]
    rcParams.update({'figure.autolayout': True})
    
    for i in metrics:
    
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        
        ax.spines['top'].set_color('w')
        ax.spines['bottom'].set_color('w')
        ax.spines['left'].set_color('w')
        ax.spines['right'].set_color('w')        
        ax.xaxis.label.set_color('w')
        ax.tick_params(axis='both', colors='w', length=8, width=4)
        
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(4)
            
        metric = i[0]
        metlist = i[1]
        
        nmacs = 0
        nmics = 0
        
        for j, k in enumerate(klist):
            plt.scatter(np.log(Nlist[j]), np.log(metlist[j]), color = k, s = 15,
            linewidths = 0.0)
            if k == 'Yellow': nmics += 1
            elif k == 'Cyan': nmacs += 1
            
        plt.scatter([0],[-1], color = 'Yellow', s=20)
        plt.scatter([0],[-1], color = 'Cyan', s=20)
        
        plt.scatter([0],[-1], color = 'Yellow', s=10,
                     label= 'microbial (n='+str(nmics)+')', alpha=0.0)
        plt.scatter([0],[-1], color = 'Cyan', s=10,
                     label= 'macrobial (n='+str(nmacs)+')', alpha=0.0)
        
        l = plt.legend(bbox_to_anchor=(0.43, .81, 0.2, .2), ncol=1,
                   prop={'size':fs-5})         
        
        l.get_frame().set_alpha(0.0)
        txtclrs = ['Yellow', 'Cyan']
        for ti, text in enumerate(l.get_texts()):
            text.set_color(txtclrs[ti])
                
        plt.xlabel('Total abundance', fontsize=fs+9, color='w')
        plt.ylabel(metric, fontsize=fs+9, color='w')
        plt.tick_params(axis='both', which='major', labelsize=fs+2)
        
        X = list(np.log(Nlist))
        Y = list(np.log(metlist))
        slope, intercept, rval, pval, stderr = sc.stats.linregress(X,Y)
        z = np.polyfit(X,Y,1)
        print 'r-squared and slope for scatterNmax:',rval**2, slope, pval
        p = np.poly1d(z)
        xp = np.linspace(min(X), max(X), 1000)
        plt.plot(xp,p(xp),'-',c='w',lw=2)
        
        plt.ylim(0, 14) #np.log(max(metlist)))
        plt.xlim(2, 14)
        
        #ax = fig.add_subplot(2, 2, 2)
    
        #for j, k in enumerate(klist):
        #    plt.scatter(np.log(Slist[j]), np.log(metlist[j]), color = k, s = 5,
        #    linewidths = 0.0)
        
        # = np.log(max(metlist))
        #plt.plot([0,xmax],[0,xmax], 'k--', lw=1)
        
        #X = list(np.log(Slist))
        #Y = list(np.log(metlist))
        #slope, intercept, rval, pval, stderr = sc.stats.linregress(X,Y)
        #z = np.polyfit(X,Y,1)
        #print 'r-squared and slope for scatterNmax:',rval**2, slope
        #p = np.poly1d(z)
        #xp = np.linspace(min(X), max(X), 1000)
        #plt.plot(xp,p(xp),'-',c='k',lw=1)
        
        #plt.xlabel('Richness', fontsize=fs)
        #plt.ylabel(metric, fontsize=fs)
        #plt.tick_params(axis='both', which='major', labelsize=fs-3)
        
        #plt.xlim(min(np.log(Slist)),max(np.log(Slist)))
        #plt.ylim(min(np.log(metlist)), max(np.log(metlist)))
        
        plt.subplots_adjust(wspace=0.4, hspace=0.4)
        
        ax.spines['top'].set_color('w')
        ax.spines['bottom'].set_color('w')
        ax.spines['left'].set_color('w')
        ax.spines['right'].set_color('w')        
        ax.xaxis.label.set_color('w')
        ax.tick_params(axis='both', colors='w')
        ax.tick_params(axis='both', colors='w', length=8, width=4)
        
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(4)
        
        rcParams.update({'figure.autolayout': True})
        plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/Nmax_scatter.png',
        transparent=True, dpi=600, pad_inches = 0.1)
        plt.close()
        
        print metric,'scatter plot: done'
    return




def Percentile():
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    
    ax.spines['top'].set_color('w')
    ax.spines['bottom'].set_color('w')
    ax.spines['left'].set_color('w')
    ax.spines['right'].set_color('w')        
    ax.xaxis.label.set_color('w')
    ax.tick_params(axis='both', colors='w', length=8, width=4)
        
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(4)
    
    
    colors = ['k', 'r', 'b', 'Cyan', 'g', 'Yellow']
    percentiles = [100, 98, 96, 94, 92, 90]
    
    for i, percentileList in enumerate(Npercs):
        
        X = list(np.log(Nlist))
        Y = list(np.log(percentileList))
        
        slope, intercept, rval, pval, stderr = sc.stats.linregress(X,Y)
        z = np.polyfit(X,Y,1)
        print 'r-squared and slope for percentile:',rval**2, slope
        
        p = np.poly1d(z)
        xp = np.linspace(min(X), max(X), 1000)
        plt.plot(xp,p(xp), '-', c=colors[i], lw=1, label=str(percentiles[i])+'th')
    
    plt.xlabel('log(total abundance)', fontsize=fs+9, color='w')
    plt.ylabel('log(abundance)', fontsize=fs+9, color='w')
    plt.tick_params(axis='both', which='major', labelsize=fs+2)
    
    #l = plt.legend(bbox_to_anchor=(-0.02, 1, 1.03, .2), loc=1, ncol=6,
    #               mode="expand", prop={'size':fs-5})
    
    l = plt.legend(loc=1, ncol=1,
                   prop={'size':fs-5})         
                            
    l.get_frame().set_alpha(0.0)
    for text in l.get_texts():
        text.set_color('w')
        
    #l = plt.legend().set_alpha(0.0)
    #for text in l.get_texts():
    #    text.set_color('w')
            
                
    #ax = fig.add_subplot(2, 2, 2)
    #for i, percentileList in enumerate(Npercs):
            
    #    X = list(np.log(Slist))
    #    Y = list(np.log(percentileList))
    #    slope, intercept, rval, pval, stderr = sc.stats.linregress(X,Y)
    #    z = np.polyfit(X,Y,1)
    #    print 'r-squared and slope for percentile:',rval**2, slope
    #    p = np.poly1d(z)
    #    xp = np.linspace(min(X), max(X), 1000)
    #    plt.plot(xp,p(xp), '-', c=colors[i], lw=1)
    
    #plt.xlabel('log(richness)', fontsize=fs)
    #plt.ylabel('log(abundance)', fontsize=fs)
    #plt.tick_params(axis='both', which='major', labelsize=fs-3)
        
    #plt.xlim(min(np.log(Slist)),max(np.log(Slist)))
    #plt.ylim(min(np.log(metlist)), max(np.log(metlist)))
    
    plt.subplots_adjust(wspace=0.3, hspace=0.3)
    
    ax.spines['top'].set_color('w')
    ax.spines['bottom'].set_color('w')
    ax.spines['left'].set_color('w')
    ax.spines['right'].set_color('w')        
    ax.xaxis.label.set_color('w')
    ax.tick_params(axis='both', colors='w', length=8, width=4)
        
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(4)    
                
    rcParams.update({'figure.autolayout': True})
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/percentile.png',
    transparent=True, dpi=600, pad_inches = 0.1)
    plt.close()
        
    print 'percentile plot: done'
    return




def rareBioScatter4(plists, Nlist, Slist, klist):
    
    fig = plt.figure()
    ct, percent = 1, 10.0
    
    for index, plist in enumerate(plists):
        ax = fig.add_subplot(2, 2, ct)
        
        nmacs = 0
        nmics = 0
        
        for j, k in enumerate(klist):
            if plist[j] > 0.0:
            
                plt.scatter(Nlist[j], plist[j]*100, color = k, s = 5,
                linewidths = 0.0)
                if k == 'Yellow': nmics += 1
                elif k == 'Cyan': nmacs += 1
        
        plt.scatter([0],[-1], color = 'Yellow', s=10,
                     label= 'microbial (n='+str(nmics)+')')
        plt.scatter([0],[-1], color = 'Cyan', s=10,
                     label= 'macrobial (n='+str(nmacs)+')')
        
        if index == 0: 
            l = plt.legend(bbox_to_anchor=(-0.03, 1.2, 2.8, .2), loc=1, ncol=2,
                               mode="expand", prop={'size':fs-15})    
            
            l.get_frame().set_alpha(0.0)
            for text in l.get_texts():
                text.set_color('w')
        
        plt.ylim(0.0, 110)
        plt.setp(ax, yticks=[20, 40, 60, 80, 100])  
        plt.xlabel('Total abundance (N)', fontsize=fs-15, color='w')
        plt.xscale('log')
        plt.ylabel('% taxa comprising\n<'+str(float(percent))+'% of N',
        fontsize=fs-15, color='w')
        plt.tick_params(axis='both', which='major', labelsize=fs-20)
        
        ax.spines['top'].set_color('w')
        ax.spines['bottom'].set_color('w')
        ax.spines['left'].set_color('w')
        ax.spines['right'].set_color('w')        
        ax.xaxis.label.set_color('w')
        ax.tick_params(axis='both', colors='w', length=8, width=4)
        
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(4)
            
        ct+= 1
        ax = fig.add_subplot(2,2, ct)
        
        for j, k in enumerate(klist):
            if plist[j] > 0.0:
                plt.scatter(Slist[j], plist[j]*100, color = k, s = 5,
                linewidths = 0.0)
                
        plt.ylim(0.0, 100)
        plt.xlim(10, 5000)
        
        plt.xlabel('Richness', fontsize=fs-15, color='w')
        plt.xscale('log')
        plt.ylabel('% taxa comprising\n<'+str(float(percent))+'% of N',
        fontsize=fs-15 , color='w')
        plt.tick_params(axis='both', which='major', labelsize=fs-15)
        
        ct+=1
        percent = percent * 0.1
        
        ax.spines['top'].set_color('w')
        ax.spines['bottom'].set_color('w')
        ax.spines['left'].set_color('w')
        ax.spines['right'].set_color('w')        
        ax.xaxis.label.set_color('w')
        ax.tick_params(axis='both', colors='w', length=8, width=4)
        
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(4)
            
    plt.subplots_adjust(wspace=0.8, hspace=0.8)    
    rcParams.update({'figure.autolayout': True})
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/RareBio_scatter.png',
    transparent=True, dpi=600, pad_inches = 0.1)
    plt.close()
    
    print 'Rare bio scatter plots: done'
    return




def rareBioScatter(plist, Nlist, Slist, klist):
    
    rcParams.update({'figure.autolayout': True})
    fig = plt.figure()
    percent = 1.0
    
    ax = fig.add_subplot(1, 1, 1)
        
    nmacs = 0
    nmics = 0
        
    for j, k in enumerate(klist):
        if plist[j] > 0.0:
            
            plt.scatter(Nlist[j], plist[j]*100, color = k, s = 20,
            linewidths = 0.0)
            if k == 'Yellow': nmics += 1
            elif k == 'Cyan': nmacs += 1
    
    plt.scatter([0],[-1], color = 'Yellow', s=10,
                 label= 'microbial (n='+str(nmics)+')', alpha=0)
    plt.scatter([0],[-1], color = 'Cyan', s=10,
                 label= 'macrobial (n='+str(nmacs)+')', alpha=0)
    
    l = plt.legend(bbox_to_anchor=(0.43, .81, 0.2, .2), ncol=1,
                   prop={'size':fs-5})    
                
    l.get_frame().set_alpha(0.0)
    txtclrs = ['Yellow', 'Cyan']
    for ti, text in enumerate(l.get_texts()):
        text.set_color(txtclrs[ti])
    
    plt.ylim(0.0, 110)
    plt.setp(ax, yticks=[20, 40, 60, 80, 100])  
    plt.xlabel('Total abundance (N)', fontsize=fs, color='w')
    plt.xscale('log')
    plt.ylabel('Rarity', fontsize=fs, color='w')
    plt.tick_params(axis='both', which='major', labelsize=fs-10)
    
    ax.spines['top'].set_color('w')
    ax.spines['bottom'].set_color('w')
    ax.spines['left'].set_color('w')
    ax.spines['right'].set_color('w')        
    ax.xaxis.label.set_color('w')
    ax.tick_params(axis='both', colors='w', length=8, width=4)
        
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(4)
    
    plt.subplots_adjust(wspace=0.8, hspace=0.8)    
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/RareBio_scatter.png',
    transparent=True, dpi=600, pad_inches = 0.1)
    plt.close()
    
    print 'Rare bio scatter plots: done'
    return




def RADfits():
    
    macrolist = ['BBS', 'NABC', 'CBC', 'GENTRY', 'MCDB']
    microlist = ['CATLIN', 'CHU', 'LAUB', 'HYDRO', 'FUNGI']
    metricNamelist = ['Geometric series', 'Zipf']
    
    fitDATA = open('/Users/lisalocey/Desktop/RareBio/R2_by_site.txt','r')

    fs = 12
    color = str()

    geomSerieslist = []
    zipfList = []
    klist = []
    
    for data in fitDATA:
        
        data_list = data.split()
        study, site, N, S, geomS, logS, zipf = data_list
        
        if S >= 10 and study in macrolist:
            
            klist.append('Cyan')
            geomSerieslist.append(float(geomS))
            zipfList.append(float(zipf))
            
        elif S >= 10 and study in microlist:
            
            klist.append('Yellow')
            geomSerieslist.append(float(geomS))
            zipfList.append(float(zipf))
            
    fitDATA.close()
    
    metlist = [geomSerieslist, zipfList]
                                                    
    fig = plt.figure()
    for i, metric in enumerate(metricNamelist):
        
        ax = fig.add_subplot(2, 2, i+1)
        ax.spines['top'].set_color('w')
        ax.spines['bottom'].set_color('w')
        ax.spines['left'].set_color('w')
        ax.spines['right'].set_color('w')        
        ax.xaxis.label.set_color('w')
        ax.tick_params(axis='both', colors='w')
        ax.tick_params(axis='both', colors='w', length=8, width=4)
        
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(4)    
        
        metricR2s = metlist[i]
        
        print metric
        
        microlist2 = []
        macrolist2 = []
        
        for ii, j in enumerate(metricR2s):
            if metricR2s[ii] >= 0.0:
                if klist[ii] == 'Yellow':
                    microlist2.append(j)
                elif klist[ii] == 'Cyan':
                    macrolist2.append(j)
                
        macrolist2 = random.sample(macrolist2, len(microlist2))
        
        n = len(microlist2)         
        D = get_kdens(microlist2)
        ymax1 = max(D[1])
        plt.plot(D[0], D[1], color = 'Yellow', lw=3,
                 label= 'microbial (n='+str(n)+')')
        
        n = len(macrolist2)
        D = get_kdens(macrolist2)
        ymax2 = max(D[1])
        plt.plot(D[0], D[1], color = 'Cyan', lw=3, 
                 label = 'macrobial (n='+str(n)+')')
        
        ymax = max(ymax1, ymax2) + 0.2*(max(ymax1, ymax2))
        plt.ylim(0.0, ymax)
        
        if i == 0:
            #l = plt.legend(bbox_to_anchor=(-0.05, 1, 2.5, .2), loc=1, ncol=2,
            #                    mode="expand", prop={'size':fs-5}) 
                
            l = plt.legend(loc=1, ncol=1,
                   prop={'size':fs-5}) 
            
            l.get_frame().set_alpha(0.0)
            for text in l.get_texts():
                text.set_color('w')
                
                
        plt.ylabel("Probability density",fontsize=fs, color='w')
        plt.xlabel(metric+' R-square', fontsize=fs, color='w')
        plt.tick_params(axis='both', which='major', labelsize=fs-3)
        
    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    rcParams.update({'figure.autolayout': True})
        
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/kdensRADfits_0trunc_noFIA.png',
    transparent=True, dpi=600, pad_inches = 0.1)
    plt.close()
    print 'RAD fits: done'     
    return



def RADfitsHist():
    
    macrolist = ['BBS', 'NABC', 'CBC', 'GENTRY', 'MCDB']
    microlist = ['CATLIN', 'CHU', 'LAUB', 'HYDRO', 'FUNGI']
    metricNamelist = ['Geometric series', 'Zipf']
    
    fitDATA = open('/Users/lisalocey/Desktop/RareBio/R2_by_site.txt', 'r')

    fs = 12
    color = str()

    geomSerieslist = []
    zipfList = []
    klist = []
    
    for data in fitDATA:
        
        data_list = data.split()
        study, site, N, S, geomS, logS, zipf = data_list
        
        if S >= 10 and study in macrolist:
            
            klist.append('Cyan')
            geomSerieslist.append(float(geomS))
            zipfList.append(float(zipf))
            
        elif S >= 10 and study in microlist:
            
            klist.append('Yellow')
            geomSerieslist.append(float(geomS))
            zipfList.append(float(zipf))
            
    fitDATA.close()
    
    metlist = [geomSerieslist, zipfList]
                                                    
    for i, metric in enumerate(metricNamelist):
        fig = plt.figure()
        
        ax = fig.add_subplot(1, 1, 1)
        metricR2s = metlist[i]
        
        print metric
        
        microlist2 = []
        macrolist2 = []
        
        for i, j in enumerate(metricR2s):
            if metricR2s[i] >= 0.0:
                if klist[i] == 'Yellow':
                    microlist2.append(j)
                elif klist[i] == 'Cyan':
                    macrolist2.append(j)
                
        macrolist2 = random.sample(macrolist2, len(microlist2))
        
        plt.hist(macrolist2,40,color='Cyan', linewidth='3', histtype='step',
        label = 'macrobial (n='+str(len(macrolist2))+')')
        
        plt.hist(microlist2,40,color='Yellow', linewidth='3', histtype='step',
        label = 'microbial (n='+str(len(microlist2))+')')
        
        l = plt.legend(bbox_to_anchor=(-0.01, 0.95, 1.02, .2), loc=1, ncol=2,
                    mode="expand", prop={'size':fs-5}) 
        
        l.get_frame().set_alpha(0.0)
        for text in l.get_texts():
            text.set_color('w')
                
        plt.xlim(0.0, 1.0)
        plt.ylabel("Number of communities",fontsize=25, color='w')
        plt.xlabel('R-square', fontsize=25, color='w')
        plt.tick_params(axis='both', which='major', labelsize=15)
        
        plt.subplots_adjust(wspace=0.4, hspace=0.4)
        
        ax.spines['top'].set_color('w')
        ax.spines['bottom'].set_color('w')
        ax.spines['left'].set_color('w')
        ax.spines['right'].set_color('w')        
        ax.xaxis.label.set_color('w')
        ax.tick_params(axis='both', colors='w', length=8, width=4)
        
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(4)
    
        rcParams.update({'figure.autolayout': True})
        plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/'+metric+'_histFits.png',
        transparent=True, dpi=600, pad_inches = 0.1)
        plt.close()
    print 'RAD fits: done'     
    return



def Fitscatter():
    
    macrolist = ['BBS', 'NABC', 'CBC', 'GENTRY', 'MCDB']
    microlist = ['CATLIN', 'CHU', 'LAUB', 'HYDRO', 'FUNGI']
    metricNamelist = ['Geometric series', 'Zipf']
    
    fitDATA = open('/Users/lisalocey/Desktop/RareBio/R2_by_site.txt','r')

    fs = 12
    color = str()

    geomSerieslist = []
    Nlist = []
    Slist = []
    zipfList = []
    klist = []
    
    for data in fitDATA:
        
        data_list = data.split()
        study, site, N, S, geomS, logS, zipf = data_list
        
        if S >= 10 and study in macrolist:
            
            klist.append('Cyan')
            geomSerieslist.append(float(geomS))
            Nlist.append(float(N))
            Slist.append(float(S))
            zipfList.append(float(zipf))
            
        elif S >= 10 and study in microlist:
            
            klist.append('Yellow')
            geomSerieslist.append(float(geomS))
            zipfList.append(float(zipf))
            Nlist.append(float(N))
            Slist.append(float(S))
            
    fitDATA.close()
    metlist = [['Geometric Series', geomSerieslist], ['Zipf', zipfList]]
                                                    
    fig = plt.figure()
    ct = 1
    for index, _list in enumerate(metlist):
        
        ax = fig.add_subplot(2, 2, ct)
        
        ax.spines['top'].set_color('w')
        ax.spines['bottom'].set_color('w')
        ax.spines['left'].set_color('w')
        ax.spines['right'].set_color('w')        
        ax.xaxis.label.set_color('w')
        ax.tick_params(axis='both', colors='w')
        ax.tick_params(axis='both', colors='w', length=8, width=4)
        
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(4)
        
        
        metricR2s = _list[1]
        name = _list[0]
        
        microlist = []
        macrolist = []
        for i, R2 in enumerate(metricR2s):
            if R2 >= 0.0:
                if klist[i] == 'Yellow':
                    microlist.append([R2, 'Yellow', Nlist[i], Slist[i]])
                    
                elif klist[i] == 'Cyan':
                    macrolist.append([R2, 'Cyan', Nlist[i], Slist[i]])
                
        macrolist = random.sample(macrolist, len(microlist))
        microlist.extend(macrolist)
        
        for j, jlist in enumerate(microlist):
            if jlist[0] > 0.0:
                plt.scatter(jlist[2], jlist[0], color = jlist[1], s = 5,
                linewidths = 0.0)
        
        """
        if index == 0: 
            l = plt.legend(bbox_to_anchor=(-0.03, 1.1, 2.47, .2), loc=1, ncol=2,
                               mode="expand", prop={'size':fs})    
            
            l.get_frame().set_alpha(0.0)
            for text in l.get_texts():
                text.set_color('w')
        """
        
        plt.ylabel(name + 'R-square',fontsize=fs, color='w')
        plt.xlabel('Total abundance', fontsize=fs, color='w')
        
        plt.ylim(0,1)
        plt.xlim(10,10**6)
        plt.xscale('log')
        plt.tick_params(axis='both', which='major', labelsize=fs-3)
        
        
        ct+=1
        ax = fig.add_subplot(2, 2, ct)
        
        ax.spines['top'].set_color('w')
        ax.spines['bottom'].set_color('w')
        ax.spines['left'].set_color('w')
        ax.spines['right'].set_color('w')        
        ax.xaxis.label.set_color('w')
        ax.tick_params(axis='both', colors='w', length=8, width=4)
        
        for axis in ['top','bottom','left','right']:
    	    ax.spines[axis].set_linewidth(4)
    
        
        for j, jlist in enumerate(microlist):
            if jlist[0] > 0.0:
                plt.scatter(jlist[3], jlist[0], color = jlist[1], s = 5,
                linewidths = 0.0)
        
        
        plt.ylabel(name + 'R-square',fontsize=fs, color='w')
        plt.xlabel('Richness', fontsize=fs, color='w')
        plt.ylim(0,1)
        plt.xlim(10,10**4)
        plt.xscale('log')
        plt.tick_params(axis='both', which='major', labelsize=fs-3)
        
        ct+=1
        
    plt.subplots_adjust(wspace=0.4, hspace=0.4)
        
    rcParams.update({'figure.autolayout': True})
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/kdensRADscatterfits_0trunc_noFIA.png',
    transparent=True, dpi=600, pad_inches = 0.1)
    plt.close()
    print 'RAD scatter fits: done'     
    return
    
        

def kdensDominance():
    metrics = [['Berger-Parker', BPlist], ['Simpsons Dominance', SimpDomlist]]

    fig = plt.figure()
    for index, i in enumerate(metrics):
    
        ax = fig.add_subplot(2, 2, index+1)

        ax.spines['top'].set_color('w')
        ax.spines['bottom'].set_color('w')
        ax.spines['left'].set_color('w')
        ax.spines['right'].set_color('w')        
        ax.xaxis.label.set_color('w')
        ax.tick_params(axis='both', colors='w')
        ax.tick_params(axis='both', colors='w', length=8, width=4)
        
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(4)
            
        
        metric = i[0]
        metlist = i[1]

        microlist = []
        macrolist = []
        
        for i, j in enumerate(metlist):
            if klist[i] == 'Yellow':
                microlist.append(j)
            elif klist[i] == 'Cyan':
                macrolist.append(j)
        
        D = get_kdens(microlist)
        n = len(microlist)
        plt.plot(D[0], D[1], color = 'Yellow', lw=3,
                 label= 'microbial (n='+str(n)+')')
        ymax1 = max(D[1])
    
        n = len(macrolist)
        D = get_kdens(macrolist)
        plt.plot(D[0], D[1], color = 'Cyan', lw=3,
                 label = 'macrobial (n='+str(n)+')')
        
        if index == 0: 
            #l = plt.legend(bbox_to_anchor=(-0.03, 1.1, 2.47, .2), loc=1, ncol=2,
            #                   mode="expand", prop={'size':fs-5})    
            
            l = plt.legend(loc=1, ncol=1,
                   prop={'size':fs-5}) 
            
            l.get_frame().set_alpha(0.0)
            for text in l.get_texts():
                text.set_color('w')

        plt.xlim(0.0, 1.0)
        plt.ylabel("Probability density",fontsize=fs, color='w')
        plt.xlabel(metric, fontsize=fs, color='w')
        plt.tick_params(axis='both', which='major', labelsize=fs)

    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    
    ax.spines['top'].set_color('w')
    ax.spines['bottom'].set_color('w')
    ax.spines['left'].set_color('w')
    ax.spines['right'].set_color('w')        
    ax.xaxis.label.set_color('w')
    ax.tick_params(axis='both', colors='w', length=8, width=4)
        
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(4)
     
    rcParams.update({'figure.autolayout': True})   
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/kdens_Dominance.png',
    transparent=True, dpi=600, pad_inches = 0.1)
    
    plt.close()
    print 'kdens fig for metrics: done'   
    
    return

    
    
def histDominance():
    metrics = [['Berger-Parker', BPlist], ['Simpsons Dominance', SimpDomlist]]
    location = 1
    fig = plt.figure()
    for index, i in enumerate(metrics):
    
        ax = fig.add_subplot(2, 2, index+1)

        metric = i[0]
        metlist = i[1]

        microlist = []
        macrolist = []
        
        for i, j in enumerate(metlist):
            if klist[i] == 'Yellow':
                microlist.append(j)
            elif klist[i] == 'Cyan':
                macrolist.append(j)
        
        plt.hist(macrolist,20,color='Cyan', linewidth='2', histtype='step',
        label = 'macrobial (n='+str(len(macrolist))+')')
        
        plt.hist(microlist,20,color='Yellow', linewidth='2', histtype='step',
        label = 'microbial (n='+str(len(microlist))+')')
        
        if index == 0:
            #l = plt.legend(bbox_to_anchor=(-0.05, 1, 2.5, .2), loc=1, ncol=2,
            #                    mode="expand", prop={'size':fs-5}) 
                
            l = plt.legend(loc=1, ncol=1,
                   prop={'size':fs-5}) 
            
            l.get_frame().set_alpha(0.0)
            for text in l.get_texts():
                text.set_color('w')
        
        location+=1
        #ymax = max(ymax1, ymax2) + 0.4*(max(ymax1, ymax2))
        #plt.ylim(0.0, ymax)
        #plt.xlim(0.0, 1.0)
        
        plt.ylabel("Number of communities",fontsize=fs, color='w')
        plt.xlabel(metric, fontsize=fs, color='w')
        plt.tick_params(axis='both', which='major', labelsize=fs)
        
        ax.spines['top'].set_color('w')
        ax.spines['bottom'].set_color('w')
        ax.spines['left'].set_color('w')
        ax.spines['right'].set_color('w')        
        ax.xaxis.label.set_color('w')
        ax.tick_params(axis='both', colors='w')            
        ax.tick_params(axis='both', colors='w', length=8, width=4)
        
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(4)
    
        
    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    
    rcParams.update({'figure.autolayout': True})    
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/hist_Dominance.png',
    transparent=True, dpi=600, pad_inches = 0.1)
    print 'hist fig for dominance: done'   
    
    return



def BPhistDominance():
    metrics = ['Berger-Parker', BPlist]
    location = 1
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    metric = metrics[0]
    metlist = metrics[1]

    microlist = []
    macrolist = []
        
    for i, j in enumerate(metlist):
        if klist[i] == 'Yellow':
            microlist.append(j)
        elif klist[i] == 'Cyan':
            macrolist.append(j)
    
    avg = np.median(macrolist)
    ser = sc.stats.sem(macrolist)
    plt.hist(macrolist,20,color='Cyan', linewidth='3', histtype='step',
    label = 'macrobial: n='+str(len(macrolist)))#+', avg BP= '+str(round(avg,2)))
    
    #print microlist
    
    avg = np.mean(microlist)
    ser = sc.stats.sem(microlist)
    label = 'microbial: n='+str(len(microlist))#+', avg BP= '+str(round(avg, 2))

    plt.hist(microlist,20,color='Yellow', linewidth='3', histtype='step', 
        label=label)
    
    #l = plt.legend(bbox_to_anchor=(-0.015, 0.95, 1.01, .2), loc=1, ncol=2,
    #                           mode="expand", prop={'size':fs-5}) 
            
    l = plt.legend(loc=1, ncol=1,
                   prop={'size':fs-5}) 
    
    l.get_frame().set_alpha(0.0)
    for text in l.get_texts():
        text.set_color('w')

    location+=1
    
    plt.xlim(0.0, 1.0)
    
    plt.ylabel("Number of communities", color='w', fontsize=fs+4)
    plt.xlabel('Dominance, '+ metric, fontsize=fs+4,  color='w')
    plt.tick_params(axis='both', which='major', labelsize=fs)

    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    
    ax.spines['top'].set_color('w')
    ax.spines['bottom'].set_color('w')
    ax.spines['left'].set_color('w')
    ax.spines['right'].set_color('w')        
    ax.xaxis.label.set_color('w')
    ax.tick_params(axis='both', colors='w', length=8, width=4)
        
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(4)
    
    rcParams.update({'figure.autolayout': True})
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/hist_Dominance.png',
    transparent=True, dpi=600, pad_inches = 0.1)
    
    #plt.show()
    #plt.close()
    print 'hist fig for Berger-Parker dominance: done'   
    
    return


def kdensSingletons():
    metrics = [['% Singletons', SinglesList]]

    fig = plt.figure()
    for index, i in enumerate(metrics):
    
        ax = fig.add_subplot(1, 1, 1)

        metric = i[0]
        metlist = i[1]

        microlist = []
        macrolist = []
        
        for i, j in enumerate(metlist):
            if klist[i] == 'Yellow':
                microlist.append(j)
            elif klist[i] == 'Cyan':
                macrolist.append(j)
        
        D = get_kdens(microlist)
        n = len(microlist)
        plt.plot(D[0], D[1], color = 'Yellow', lw=6,
                 label= 'microbial (n='+str(n)+')')
        ymax1 = max(D[1])
    
        n = len(macrolist)
        D = get_kdens(macrolist)
        plt.plot(D[0], D[1], color = 'Cyan', lw=6,
                 label = 'macrobial (n='+str(n)+')')
        ymax2 = max(D[1])
    
        if index == 0: 
            #l = plt.legend(bbox_to_anchor=(-0.03, 1.1, 2.47, .2), loc=1, ncol=2,
            #                   mode="expand", prop={'size':fs-5})    
            
            l = plt.legend(loc=1, ncol=1,
                   prop={'size':fs-5}) 
            
            l.get_frame().set_alpha(0.0)
            for text in l.get_texts():
                text.set_color('w')
        
        ymax = max(ymax1, ymax2) + 0.4*(max(ymax1, ymax2))
        plt.ylim(0.0, ymax)
        plt.xlim(0.0, 1.0)
        plt.ylabel("Probability density",fontsize=fs+3, color='w')
        plt.xlabel(metric, fontsize=fs+3, color='w')
        plt.tick_params(axis='both', which='major', labelsize=fs)

    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    
    ax.spines['top'].set_color('w')
    ax.spines['bottom'].set_color('w')
    ax.spines['left'].set_color('w')
    ax.spines['right'].set_color('w')        
    ax.xaxis.label.set_color('w')
    ax.tick_params(axis='both', colors='w', length=8, width=4)
        
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(4)
    
    rcParams.update({'figure.autolayout': True})    
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/kdens_Singletons.png',
    transparent=True, dpi=600, pad_inches = 0.1)
    
    plt.close()
    print 'kdens fig for singletons: done'   


def rarebioScale(plists, klist):
    
    rcParams.update({'figure.autolayout': True})
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    
    Xs = [0.01, 0.1, 1.0, 10.0, 80.0]
    
    microMeans = []
    microLow = []
    microHigh = []
    microSEs = []
    
    macroMeans = []
    macroLow = []
    macroHigh = []
    macroSEs = []
    
    for j, plist in enumerate(plists):
        
        macrolist = []
        microlist = []
        
        for i, k in enumerate(plist):
            if klist[i] == 'Yellow':
                microlist.append(k*100)
            elif klist[i] == 'Cyan':
                macrolist.append(k*100)
        
        macroMeans.append(np.mean(macrolist))
        macroLow.append(np.percentile(macrolist, 10))
        macroHigh.append(np.percentile(macrolist, 90))
        macroSEs.append(sem(macrolist, ddof=1))
        
        microMeans.append(np.mean(microlist))
        microLow.append(np.percentile(microlist, 10))
        microHigh.append(np.percentile(microlist, 90))
        microSEs.append(sem(microlist, ddof=1))
        
    plt.fill_between(Xs, macroLow, macroHigh, color='Cyan', linewidth=0.0, label='macrobes', alpha=0.6) 
    plt.plot(Xs, macroMeans, color='Cyan', ls='--', lw=3, label='macrobes')
    #plt.errorbar(Xs, macroMeans, yerr=macroSEs, lw=3, color = 'Cyan')
    
    plt.fill_between(Xs, microLow, microHigh, color='Yellow', linewidth=0.0, label='microbes', alpha=0.6) 
    plt.plot(Xs, microMeans, color='Yellow', ls='--', label='microbes', lw=3)
    #plt.errorbar(Xs, microMeans, yerr=microSEs, lw=3, color = 'Yellow')
    
            
    l = plt.legend(loc=2, ncol=1,
                   prop={'size':fs-5}) 
    
    l.get_frame().set_alpha(0.0)
    for text in l.get_texts():
        text.set_color('w')
    
    
    plt.ylim(0,100)
    plt.xscale('log')
    plt.ylabel("% taxa",fontsize=fs, color='w')
    plt.xlabel('% Total abundance', fontsize=fs, color='w')
    
    ax.spines['top'].set_color('w')
    ax.spines['bottom'].set_color('w')
    ax.spines['left'].set_color('w')
    ax.spines['right'].set_color('w')        
    ax.xaxis.label.set_color('w')
    ax.tick_params(axis='both', colors='w', length=8, width=4)
        
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(4)
    
    plt.tick_params(axis='both', which='major', labelsize=fs-10)
    plt.subplots_adjust(wspace=0.8, hspace=0.8)
    
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/rarebioScale_hulls.png',
    transparent=True, dpi=600, pad_inches = 0.3)
    
    plt.close()
    
    
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    plt.errorbar(Xs, microMeans, yerr=microSEs, lw=3, color='Yellow')
    plt.errorbar(Xs, macroMeans, yerr=macroSEs, lw=3, color='Cyan')
    
    plt.scatter([0],[-1], color = 'Yellow',
                 label= 'microbes', alpha=0)
    plt.scatter([0],[-1], color = 'Cyan',
                 label= 'macrobes', alpha=0)
    
    l = plt.legend(loc=2, ncol=1,
                   prop={'size':fs-5}) 
    
    l.get_frame().set_alpha(0.0)
    txtclrs = ['Yellow', 'Cyan']
    for ti, text in enumerate(l.get_texts()):
        text.set_color(txtclrs[ti])
    
    plt.ylim(0,100)
    plt.xscale('log')    
    plt.ylabel("% taxa",fontsize=fs, color='w')
    plt.xlabel('% Total abundance', fontsize=fs, color='w')
    
    ax.spines['top'].set_color('w')
    ax.spines['bottom'].set_color('w')
    ax.spines['left'].set_color('w')
    ax.spines['right'].set_color('w')        
    ax.xaxis.label.set_color('w')
    ax.tick_params(axis='both', colors='w', length=8, width=4)
        
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(4)
    
    
    plt.tick_params(axis='both', which='major', labelsize=fs-10)
    plt.subplots_adjust(wspace=0.8, hspace=0.8)
    
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/rarebioScale_error.png',
    transparent=True, dpi=600, pad_inches = 0.3)
    
    print 'kdens for rarebioScale: done'
    plt.close()
    
    return

    
    
    
def domTailScale(dlists, klist):
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    Xs = [1, 2, 3, 4, 5]
    
    microMeans = []
    microLow = []
    microHigh = []
    microSEs = []
    
    macroMeans = []
    macroLow = []
    macroHigh = []
    macroSEs = []
    
    for j, dlist in enumerate(dlists):
        
        macrolist = []
        microlist = []
        
        for i, k in enumerate(dlist):
            if klist[i] == 'Yellow':
                microlist.append(k*100)
            elif klist[i] == 'Cyan':
                macrolist.append(k*100)
        
        macroMeans.append(np.mean(macrolist))
        macroLow.append(np.percentile(macrolist, 10))
        macroHigh.append(np.percentile(macrolist, 90))
        macroSEs.append(sem(macrolist, ddof=1))
        
        microMeans.append(np.mean(microlist))
        microLow.append(np.percentile(microlist, 10))
        microHigh.append(np.percentile(microlist, 90))
        microSEs.append(sem(microlist, ddof=1))
    
        
    plt.fill_between(Xs, macroLow, macroHigh, color='Cyan', linewidth=0.0) 
    plt.plot(Xs, macroMeans, color='Cyan', ls='--', lw=2)
    #plt.errorbar(Xs, macroMeans, yerr=macroSEs, lw=3, color = 'Cyan')
    
    plt.fill_between(Xs, microLow, microHigh, color='Yellow', linewidth=0.0) 
    plt.plot(Xs, microMeans, color='Yellow', ls='--', lw=2)
    #plt.errorbar(Xs, microMeans, yerr=microSEs, lw=3, color = 'Yellow')
    
    plt.tick_params(axis='both', which='major', labelsize=fs)
    plt.setp(ax, xticks=[1, 2, 3, 4, 5])
    plt.xlim(1,5)
    plt.ylim(0,100)
    
    plt.ylabel("% Total abundance",fontsize=fs+6, color='w')
    plt.xlabel('no. dominant taxa', fontsize=fs+6, color='w')
    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    
    ax.spines['top'].set_color('w')
    ax.spines['bottom'].set_color('w')
    ax.spines['left'].set_color('w')
    ax.spines['right'].set_color('w')        
    ax.xaxis.label.set_color('w')
    ax.tick_params(axis='both', colors='w', length=8, width=4)
        
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(4)
    
    rcParams.update({'figure.autolayout': True})
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/domTailScale_hulls.png',
    transparent=True, dpi=600, pad_inches = 0.1)
        
    plt.close()
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.errorbar(Xs, macroMeans, yerr=macroSEs, lw=3, color = 'Cyan', label='macrobes')
    plt.errorbar(Xs, microMeans, yerr=microSEs, lw=3, color = 'Yellow', label='microbes')
    
    #l = plt.legend(bbox_to_anchor=(-0.01, 0.95, 1, .2), loc=1, ncol=2,
    #               mode="expand", prop={'size':fs-5})    
            
    l = plt.legend(loc=1, ncol=1,
                   prop={'size':fs-5}) 
                   
    l.get_frame().set_alpha(0.0)
    for text in l.get_texts():
        text.set_color('w')     
    
    plt.ylabel("% Total abundance",fontsize=fs+6, color='w')
    plt.xlabel('no. dominant taxa', fontsize=fs+6, color='w')
    plt.tick_params(axis='both', which='major', labelsize=fs)
    plt.setp(ax, xticks=[1, 2, 3, 4, 5])
    plt.xlim(1,5)
    plt.ylim(0,100)
    
    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    
    ax.spines['top'].set_color('w')
    ax.spines['bottom'].set_color('w')
    ax.spines['left'].set_color('w')
    ax.spines['right'].set_color('w')        
    ax.xaxis.label.set_color('w')
    ax.tick_params(axis='both', colors='w', length=8, width=4)
        
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(4)
        
    rcParams.update({'figure.autolayout': True})
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/domTailScale_error.png',
    transparent=True, dpi=600, pad_inches = 0.1)
    
    print 'kdens for dominant taxa: done'
    plt.close()
    
    return




def NandSkdens(Nlist, Slist, klist):

    metrics = [['log(N)', Nlist], ['log(S)', Slist]]

    fig = plt.figure()
    for index, i in enumerate(metrics):
    
        ax = fig.add_subplot(2,2,index+1)
        
        ax.spines['top'].set_color('w')
        ax.spines['bottom'].set_color('w')
        ax.spines['left'].set_color('w')
        ax.spines['right'].set_color('w')        
        ax.xaxis.label.set_color('w')
        ax.tick_params(axis='both', colors='w', length=8, width=4)
        
        for axis in ['top','bottom','left','right']:
    	    ax.spines[axis].set_linewidth(4)
    
        metric = i[0]
        metlist = np.log(i[1])
        
        microlist = []
        macrolist = []
        
        for i, j in enumerate(metlist):
            if klist[i] == 'Yellow':
                microlist.append(j)
            elif klist[i] == 'Cyan':
                macrolist.append(j)
        
        D = get_kdens(microlist)
        n = len(microlist)
        plt.plot(D[0], D[1], color = 'Yellow', lw=3,
                 label= 'microbial (n='+str(n)+')')
        ymax1 = max(D[1])
    
        n = len(macrolist)
        D = get_kdens(macrolist)
        plt.plot(D[0], D[1], color = 'Cyan', lw=3,
                 label = 'macrobial (n='+str(n)+')')
        ymax2 = max(D[1])
    
        if index == 0: 
            l = plt.legend(bbox_to_anchor=(-0.03, 1.1, 2.47, .2), loc=1, ncol=2,
                               mode="expand", prop={'size':fs-15})    
            
            l.get_frame().set_alpha(0.0)
            for text in l.get_texts():
                text.set_color('w')

        
        #ymax = max(ymax1, ymax2) + 0.4*(max(ymax1, ymax2))
        #plt.ylim(0.0, ymax)
        #plt.xlim(0.0, 1.0)
        plt.ylabel("Probability density",fontsize=fs-15, color='w')
        plt.xlabel(metric, fontsize=fs-15, color='w')
        plt.tick_params(axis='both', which='major', labelsize=fs-15)

    plt.subplots_adjust(wspace=0.6, hspace=0.6)
    
    ax.spines['top'].set_color('w')
    ax.spines['bottom'].set_color('w')
    ax.spines['left'].set_color('w')
    ax.spines['right'].set_color('w')        
    ax.xaxis.label.set_color('w')
    ax.tick_params(axis='both', colors='w', length=8, width=4)
        
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(4)
        
    rcParams.update({'figure.autolayout': True})
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/NandS_kdens.png',
    transparent=True, dpi=600, pad_inches = 0.1)
    plt.close()
    print 'kdens fig for N and S: done'   



Percentile()
kdensMetrics()
histMetrics()

scatterMetrics()
scatterDominance()
kdensDominance()

BPhistDominance()
histDominance()
scatterNmax()

kdensSingletons() # not fixed
RADfits()
RADfitsHist()
Fitscatter()
    
plists = [ptzonelist, ptonelist, onelist, tenlist, atylist]    
rarebioScale(plists, klist)

dlists = [BP1list, BP2list, BP3list, BP4list, BP5list]
domTailScale(dlists, klist)

NandSkdens(Nlist, Slist, klist)
rareBioScatter(onelist, Nlist, Slist, klist)