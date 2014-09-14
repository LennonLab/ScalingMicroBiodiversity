from __future__ import division
import  matplotlib.pyplot as plt
import numpy as np
import sys
import os
import feasible_functions as ff


radDATA = open('/Users/lisalocey/Desktop/XRADs/RADdata.txt','r')

fs = 12
color = str()

Nlist = []
Slist = []
Evarlist = [] 
ESimplist = []
ENeelist = []
EHeiplist = []
EQlist = []
klist = []

onelist = []
ptonelist = []
ptzonelist = []

ct = 0
klist = []
for data in radDATA:
    break
    data_list = data.split()
    name, kind, N, S, Evar, ESimp, ENee, EHeip, EQ, one, ptone, ptzone = data_list
    
    if float(ESimp) > 1.0 or float(ESimp) < 0.0:
        print 'Simpsons evevnness violation',ESimp, name
        sys.exit()   
    
    Nlist.append(float(N))
    Slist.append(float(S))
    Evarlist.append(float(Evar))
    ESimplist.append(float(ESimp))
    ENeelist.append(float(ENee))
    EHeiplist.append(float(EHeip))
    EQlist.append(float(EQ))
    onelist.append(float(one))
    ptonelist.append(float(ptone))
    ptzonelist.append(float(ptzone))

    if kind == 'micro':
        color = 'r'
    elif kind == 'macro':
        color = 'k'
    klist.append(color)
    
    ct+=1
    #if ct >= 5000: break
    print 'community',ct,name
    
#sys.exit()    
#print 'There are',len(Nlist),'sites of ecological communities'

metrics = [['Evar', Evarlist], ['Simpsons', ESimplist], ['Eq', EQlist]]
metrics = [['Simpsons', ESimplist]]
metrics = [['Slope of RAD', ENeelist], ['Eq', EQlist]]


ct = 1
for i in metrics:
    
    break         
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    metric = i[0]
    metlist = i[1]

    metlist2, klist2 = [list(x) for x in zip(*sorted(zip(metlist, klist),
                    key=lambda pair: pair[0]))]

    metlist2.reverse()
    klist2.reverse()

    for i, val in enumerate(metlist2):
    
        size = float()
        a = float()
        if klist2[i] == 'r':
            size = 20
            a = 0.9
        else: 
            size = 5
            a = 0.9

        plt.scatter([i+1], [val], color = klist2[i], s = size, linewidths = 0.0, alpha = a)

    plt.xlabel('Rank in'+metric, fontsize=fs)
    plt.ylabel(metric, fontsize=fs)
    plt.tick_params(axis='both', which='major', labelsize=fs-3)
    plt.xlim(0,18000)
    plt.savefig(metric+'Ranks.png', dpi=600)

    print 'fig'+str(ct)+': done'   
    ct += 1


    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for j, k in enumerate(klist):
    
        size = float()
        a = float()
        if klist[i] == 'r':
            size = 15
            a = 0.6
        else: 
            size = 15
            a = 0.2
    
        plt.scatter(Nlist[j], metlist[j], color = k, alpha = a, s = size, linewidths = 0.0)

    plt.xlabel('Total abundance', fontsize=fs)
    plt.ylabel(metric, fontsize=fs)
    plt.tick_params(axis='both', which='major', labelsize=fs-3)
    plt.xscale('log')
    plt.savefig('Nv'+metric+'.png', dpi=600)

    print 'fig'+str(ct)+': done'
    ct += 1


    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for j, k in enumerate(klist):
    
        size = float()
        a = float()
        if klist[i] == 'r':
            size = 15
            a = 0.6
        else: 
            size = 15
            a = 0.2
    
        plt.scatter(Slist[j], metlist[j], color = k, alpha = a, s = size, linewidths = 0.0)

    plt.xlabel('Richness', fontsize=fs)
    plt.xscale('log')
    plt.ylabel(metric, fontsize=fs)
    plt.tick_params(axis='both', which='major', labelsize=fs-3)
    plt.savefig('Sv'+metric+'.png', dpi=600)

    print 'fig'+str(ct)+': done'
    ct += 1
    
    
    

def rarebio(plists, klist):

    fig = plt.figure()
    plists = [onelist, ptonelist, ptzonelist]    
            
    percent = 1.0
    for j, plist in enumerate(plists):
    
        ax = fig.add_subplot(2,2,j+1)
        plist2, klist2 = [list(x) for x in zip(*sorted(zip(plist, klist),
                    key=lambda pair: pair[0]))]

        plist2.reverse()
        klist2.reverse()

        for i, val in enumerate(plist2):
    
            size = float()
            a = float()
            if klist2[i] == 'r':
                size = 20
                a = 0.9
            else: 
                size = 10
                a = 0.9
        
            plt.scatter([i+1], [val], color = klist2[i], s = size, linewidths = 0.0, alpha = a)
    
        plt.xscale('log')
        plt.xlabel('log(rank)', fontsize=fs)
        plt.ylabel('Portion taxa comprising\n<'+str(percent)+' abundance', fontsize=fs)
        plt.tick_params(axis='both', which='major', labelsize=fs-3)
        #plt.xlim(0,18000)
        percent = 0.1*percent

    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    plt.savefig('/Users/lisalocey/Desktop/XRADs/RareBio.png', dpi=600)
    #plt.show()
    return


def RADfits():
    
    macrolist = ['BBS', 'FIA', 'NABC', 'CBC', 'GENTRY', 'MCDB']
    microlist = ['CATLIN', 'CHU', 'LAUB', 'HYDRO', 'FUNGI']
    metricNamelist = ['Geometric series', 'Log-series', 'Zipf']
    
    fitDATA = open('/Users/lisalocey/Desktop/XRADs/R2_by_site.txt','r')

    fs = 12
    color = str()

    logSerieslist = []
    geomSerieslist = []
    zipfList = []
    klist = []
    
    for data in fitDATA:
        
        data_list = data.split()
        study, site, N, S, geomS, logS, zipf = data_list
        
        if S >= 10 and study in microlist:
            
            klist.append('r')
            logSerieslist.append(float(logS))
            geomSerieslist.append(float(geomS))
            zipfList.append(float(zipf))
            
        elif S >= 10:
            
            klist.append('k')
            logSerieslist.append(float(logS))
            geomSerieslist.append(float(geomS))
            zipfList.append(float(zipf))
        
            
    fitDATA.close()
    
    metlist = [geomSerieslist, logSerieslist, zipfList]
                                                    
    fig = plt.figure()
    for i, metric in enumerate(metricNamelist):
        
        ax = fig.add_subplot(2,2,i+1)
        metricR2s = metlist[i]
        
        print metric
        
        metricR2s2, klist2 = [list(x) for x in zip(*sorted(zip(metricR2s, klist),
                    key=lambda pair: pair[0]))]

        metricR2s2.reverse()
        klist2.reverse()

        for i, val in enumerate(metricR2s2):
            
            print i+1
            
            size = float()
            a = float()
            if klist2[i] == 'r':
                size = 20
                a = 0.9
            else: 
                size = 5
                a = 0.9
            
            #if i >= 5000: break
            
            plt.scatter([i+1], [val], color = klist2[i], s = size, linewidths = 0.0, alpha = a)

        plt.xlabel('Rank in R-square', fontsize=fs)
        plt.ylabel(metric+' R-square', fontsize=fs)
        plt.tick_params(axis='both', which='major', labelsize=fs-3)
        #plt.ylim(0,1)

    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    plt.savefig('/Users/lisalocey/Desktop/XRADs/figs/RADfits_notrunc.png', dpi=600)
    
    print 'RAD fits: done'     
    return

RADfits()
#rarebio(plists, klist)