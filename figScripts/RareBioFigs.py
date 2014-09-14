from __future__ import division
import  matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
import random
import scipy as sc

import sys

sys.path.append("/Users/lisalocey/Desktop/RareBio/")
import ModelsN

sys.path.append("/Users/lisalocey/Desktop/RareBio/")
import feasible_functions as ff

sys.path.append("/Users/lisalocey/Desktop/RareBio/global/GenAnalysis/tools/")
import mete
#import pln

sys.path.append("/Users/lisalocey/Desktop/partitions/partitions/Python/")
import partitions



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




def RandCompFast(Q, N, sample_size, maxn=None):
    
    comps = []
    while len(comps) < sample_size:

        composition = []    
        indices = []
        
        while len(indices) < N-1:
            index = random.randint(1, Q-1)
            if index in indices: continue
            else: indices.append(index)
    
        indices.sort()
        indices.append(Q)
    
        nsum = 0
        composition = [] 
        for i in indices:
            i -= sum(composition) 
            composition.append(i)
            nsum += i
        
        composition.sort()
    	composition.reverse()
        comps.append(composition)

    return comps
    
    

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
            
            

def kdensMetrics():
    
    metrics = [['Evar', Evarlist], ['Simpson\'s Diversity, 1-D', SimpDomlist], ['Simpsons\' evenness', ESimplist], ['Berger-Parker', BPlist]]
    metrics = [['Heip\'s', EHeiplist], ['Simpson\'s Diversity, 1-D', SimpDomlist], ['Simpsons\' evenness', ESimplist], ['Berger-Parker', BPlist]]

    fig = plt.figure()
    for index, i in enumerate(metrics):
    
        fig.add_subplot(2,2,index+1)

        metric = i[0]
        metlist = i[1]

        microlist = []
        macrolist = []
        
        for i, j in enumerate(metlist):
            if klist[i] == 'm':
                microlist.append(j)
            elif klist[i] == 'c':
                macrolist.append(j)
            
        print len(microlist), len(macrolist)
                    
        D = get_kdens(microlist)
        n = len(microlist)
        plt.plot(D[0], D[1], color = 'm', lw=3, alpha = 0.9,
                 label= 'microbial (n='+str(n)+')')
        ymax1 = max(D[1])
    
        
        n = len(macrolist)
        D = get_kdens(macrolist)
        plt.plot(D[0], D[1], color = 'c', lw=3, alpha = 0.9,
                 label = 'macrobial (n='+str(n)+')')
        ymax2 = max(D[1])
    
        
        if index == 0: 
            plt.legend(bbox_to_anchor=(-0.03, 1.1, 2.47, .2), loc=10, ncol=2,
                               mode="expand",prop={'size':fs-2})    
        
        
        ymax = max(ymax1, ymax2) + 0.4*(max(ymax1, ymax2))
        plt.ylim(0.0, ymax)
        plt.ylabel("Probability density",fontsize=fs)
        plt.xlabel(metric, fontsize=fs)
        plt.tick_params(axis='both', which='major', labelsize=fs-3)

    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/ProjectN/MetricsKdensEMP.png',
    dpi=600)
    plt.close()
    print 'kdens fig for metrics: done'   
    
    return
    
    
    
    
def scatterMetrics():
    #metrics = [['Evar', Evarlist], ['Simpsons', ESimplist]]
    #metrics = [['Evar', Evarlist], ['Simpsons', ESimplist], ['Shannons', Shanlist], ['Eq', EQlist]]
    #metrics = [['Evenness (Evar)', Evarlist], ['Simpson\'s Diversity', SimpDomlist], ['Simpsons\' evenness', ESimplist], ['Berger-Parker', BPlist]]
    #metrics = [['Evar', Evarlist], ['Simpson\'s Diversity, 1-D', SimpDomlist], ['Simpsons\' evenness', ESimplist], ['Berger-Parker', BPlist]]
    metrics = [['Heip\'s', EHeiplist], ['Simpson\'s Diversity, 1-D', SimpDomlist], ['Simpsons\' evenness', ESimplist], ['Berger-Parker', BPlist]]
    
    fig = plt.figure()
    for index, i in enumerate(metrics):
    
        fig.add_subplot(2, 2, index+1)

        metric = i[0]
        metlist = i[1]
        
        MacListX = []
        MacListY = []
        
        MicListX = []
        MicListY = []
        
        nmacs = 0
        nmics = 0
        
        for j, k in enumerate(klist):
            
            if k == 'm': 
                nmics += 1
                MicListX.append(Nlist[j])
                MicListY.append(metlist[j])
                
            elif k == 'c':
                nmacs += 1
                MacListX.append(Nlist[j])
                MacListY.append(metlist[j])
                    
        plt.scatter(MacListX, MacListY, color = 'c', alpha=0.6 , s=0.1)                                                
        plt.scatter(MicListX, MicListY, color = 'm', alpha=0.6 , s=0.2)
                                                                                                                                                
        if index+1 == 1:
            plt.scatter([0],[-1], color = 'm', alpha = 0.8, s=10,
                        label= 'microbial (n='+str(nmics)+')')
            
            plt.scatter([0],[-1], color = 'c',alpha=0.8, s=10,
                        label= 'macrobial (n='+str(nmacs)+')')
            
            plt.legend(bbox_to_anchor=(-0.03, 1.1, 2.47, .2), loc=10, ncol=2,
                                mode="expand",prop={'size':fs-2})
        
        
        plt.ylim(0, 1)     
        plt.xlabel('Total abundance', fontsize=fs)
        plt.ylabel(metric, fontsize=fs)
        plt.tick_params(axis='both', which='major', labelsize=fs-3)
        plt.xscale('log')
        
        print metric,'scatter plot: done'
    
    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/ProjectN/ScatterMetrics.png', dpi=600)
    plt.close()
    
    return



def scatterMetricsRare():
    global rareSumOnesList
    List = np.array(rareSumOnesList) * 100
    rareSumOnesList = List.tolist()
    
    metrics = [['% abundance of singletons', rareSumOnesList], ['% abundance of\nrarest taxon', rareRelList], ['number of singletons', rareOnesList], ['Chance that the next sampled\nindividual is from a taxon\n of lesser abundance', rarePairList]]
     
    
    fig = plt.figure()
    for index, i in enumerate(metrics):
    
        fig.add_subplot(2, 2, index+1)

        metric = i[0]
        metlist = i[1]
        
        MacListX = []
        MacListY = []
        
        MicListX = []
        MicListY = []
        
        nmacs = 0
        nmics = 0
        
        for j, k in enumerate(klist):
            
            if k == 'm': 
                nmics += 1
                MicListX.append(Nlist[j])
                MicListY.append(metlist[j])
                
            elif k == 'c':
                nmacs += 1
                MacListX.append(Nlist[j])
                MacListY.append(metlist[j])
                
            
        plt.scatter(MacListX, MacListY, color = 'c', alpha=0.6 , s=0.1)                                                
        plt.scatter(MicListX, MicListY, color = 'm', alpha=0.6 , s=0.2)
                                                                                                                                                
        if index+1 == 1:
            plt.scatter([0],[-1], color = 'm', alpha = 0.8, s=10,
                        label= 'microbial (n='+str(nmics)+')')
            
            plt.scatter([0],[-1], color = 'c',alpha=0.8, s=10,
                        label= 'macrobial (n='+str(nmacs)+')')
            
            plt.legend(bbox_to_anchor=(-0.03, 1.1, 2.47, .2), loc=10, ncol=2,
                                mode="expand",prop={'size':fs-2})
        
        
        if index == 0:
            plt.ylim(.001, 100.0)
            plt.yscale('log')
        if index == 1:
            plt.ylim(.000001, 10.0)
            plt.yscale('log')
        if index == 2:
            plt.ylim(1, 10**5)
            plt.yscale('log')
        if index == 3:
            plt.ylim(0.0, 0.6)     
        
        
        plt.xlabel('Total abundance', fontsize=fs)
        plt.ylabel(metric, fontsize=fs)
        plt.tick_params(axis='both', which='major', labelsize=fs-3)
        plt.xscale('log')
        
        print metric,'scatter plot: done'
    
    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/ProjectN/ScatterRare.png',
    dpi=600)
    plt.close()
    
    return
    




def scatterNmax():
    
    metrics = [['Most abundant species', NmaxList]]
    
    for i in metrics:
    
        fig = plt.figure()
        fig.add_subplot(2, 2, 1)

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
        
        fig.add_subplot(2, 2, 2)
    
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



def scatterNmax2():
    
    metrics = [['Most abundant species', NmaxList]]
    
    for i in metrics:
    
        fig = plt.figure()
        

        metric = i[0]
        metlist = i[1]
        
        MacListX = []
        MacListY = []
        
        MicListX = []
        MicListY = []
        
        nmacs = 0
        nmics = 0
        i = 0
        for j, k in enumerate(klist):
            
            if k == 'm': 
                nmics += 1
                MicListX.append(Nlist[j])
                MicListY.append(metlist[j])
                i = 1
                
            elif k == 'c':
                nmacs += 1
                MacListX.append(Nlist[j])
                MacListY.append(metlist[j])
                i = 2
            
                
        fig.add_subplot(2, 2, 1) 
        
        X = list(np.log(MicListX))
        Y = list(np.log(MicListY))
        slope, intercept, rval, pval, stderr = sc.stats.linregress(X,Y)
        slope, r2 = round(slope,2), round(rval**2,2)
        print 'r-squared and slope for scatterNmax Micro:',rval**2, slope
        
        plt.plot([1,max(MicListX)],[1,max(MicListX)],c='k',lw=0.5)
        
        plt.scatter(MicListX, MicListY, color = 'm', alpha = 0.6, s = 0.5,
                     label= 'microbial (n='+str(nmics)+')\nR2='+str(r2)+' slope='+str(slope))
        
        leg = plt.legend(loc=2,prop={'size':fs-2})
        leg.draw_frame(False)

        #plt.xlim(1,max(MicListX))
        #plt.ylim(1,max(MicListX))
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Sample abundance', fontsize=fs)
        plt.ylabel(metric, fontsize=fs)
        plt.tick_params(axis='both', which='major', labelsize=fs-3)
        
        
        
        fig.add_subplot(2, 2, 2)                                                                                                                     
        X = list(np.log(MacListX))
        Y = list(np.log(MacListY))
        slope, intercept, rval, pval, stderr = sc.stats.linregress(X,Y)
        slope, r2 = round(slope,2), round(rval**2,2)
        
        plt.plot([1,max(MacListX)],[1,max(MacListX)],c='k',lw=0.5)
        print 'r-squared and slope for scatterNmax macro:',rval**2, slope
        
        plt.scatter(MacListX, MacListY, color = 'c', alpha=0.5 , s=0.5,
                     label= 'macrobial (n='+str(nmacs)+')\nR2='+str(r2)+' slope='+str(slope))
        
        leg = plt.legend(loc=2,prop={'size':fs-2})
        leg.draw_frame(False)
        
        plt.xlim(1, max(MacListX))
        plt.ylim(1, max(MacListX))
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Sample abundance', fontsize=fs)
        plt.ylabel(metric, fontsize=fs)
        plt.tick_params(axis='both', which='major', labelsize=fs-3)
        
                
    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/ScatterMax2.png',
    dpi=600)
    plt.close()
        
    return





def Fig1(tool):
    
    fig = plt.figure()
    
    klist1 = []
    Nlist1 = []
    Slist1 = []
    ObsHeip1 = []
    ObsBP1 = []
    ObsRare1 = []
    
    for i, j in enumerate(klist):
        if j == 'm' or j == 'c': 
            
            klist1.append(j)
            Nlist1.append(Nlist[i])
            Slist1.append(Slist[i])
            ObsHeip1.append(EHeiplist[i])
            ObsBP1.append(BPlist[i])
            ObsRare1.append(rareOnesList[i])
            
    Indices = range(len(Nlist1))
    
    #IndicesSample = random.sample(Indices, 200) # w/out replacement
    IndicesSample = list(Indices)
    #random.shuffle(IndicesSample)
    
    NlistSample1 = []
    klistSample1 = []
    SlistSample1 = []
    ObsHeipS1 = []
    ObsBPS1 = []
    ObsRareS1 = []
    
    
    for i in IndicesSample:
        
        #if len(NlistSample1) > 50: break
        if Nlist1[i] < 10**4: continue
        if Nlist1[i] > 10**5: continue
        if Slist[i] > 300: continue
        
        NlistSample1.append(Nlist1[i])
        klistSample1.append(klist1[i])
        SlistSample1.append(Slist1[i])
        ObsHeipS1.append(ObsHeip1[i])
        ObsBPS1.append(ObsBP1[i])
        ObsRareS1.append(ObsRare1[i])
    
    FS_RADs = ff.getFS(NlistSample1, SlistSample1, tool)    
    sys.exit()
                            
    IN = open('/Users/lisalocey/Desktop/RareBio/macrostates/'+tool+'.txt','r')
    
    RADdict = {}
    for line in IN:
        
        RAD = eval(line)
        n = sum(RAD)
        s = len(RAD)
        nsLabel = (n, s)
        if nsLabel in RADdict: RADdict[nsLabel].append(RAD)
        else: RADdict[nsLabel] = [RAD]
    
    IN.close()
    
    
    NlistSample = []
    klistSample = []
    SlistSample = []
    
    ObsHeip = []
    ObsBP = []
    ObsRare = []
    
    FS_Heip = []
    FS_BP = []
    FS_Rare = []
    
    for i, N in enumerate(NlistSample1):
        
        N = int(N)
        S = int(SlistSample1[i])
        
        
        nsLabel = (N, S)
        
        if nsLabel in RADdict:
            NlistSample.append(N)
            SlistSample.append(S)
            klistSample.append(klistSample1[i])
            
            ObsHeip.append(ObsHeipS1[i])
            ObsBP.append(ObsBPS1[i])
            ObsRare.append(ObsRareS1[i])
            
            FS_RADs = RADdict[(N, S)]
            FS_HeipList = []
            FS_BPList = []
            FS_RareList = []
    
            for RAD in FS_RADs:
            
                Heip = ff.Heips_evenness(RAD)
                FS_HeipList.append(Heip)
        
                BP = ff.Berger_Parker(RAD)
                FS_BPList.append(BP)
        
                Rare = ff.rarityOnes(RAD)  
                FS_RareList.append(Rare)
              
            FS_Heip.append(np.mean(FS_HeipList))
            FS_BP.append(np.mean(FS_BPList))
    	    FS_Rare.append(np.mean(FS_RareList))
    
    print len(EHeiplist), len(FS_Heip)
    
    
    fig.add_subplot(3, 3, 1)
    """ kdens of Heips """
    
    metric = 'Evenness, Heip\'s'
    microlist = []
    macrolist = []
        
    for i, j in enumerate(EHeiplist):
        if klist[i] == 'm':
            microlist.append(j)
        elif klist[i] == 'c':
            macrolist.append(j)
            
    D = get_kdens(microlist)
    n = len(microlist)
    plt.plot(D[0], D[1], color = 'm', lw=3, alpha = 0.9,
             label= 'microbial (n='+str(n)+')')
    ymax1 = max(D[1])
    
        
    n = len(macrolist)
    D = get_kdens(macrolist)
    plt.plot(D[0], D[1], color = 'c', lw=3, alpha = 0.9,
             label = 'macrobial (n='+str(n)+')')
    
    ymax2 = max(D[1])
    ymax = max(ymax1, ymax2) + 0.4*(max(ymax1, ymax2))
    
    plt.ylim(0.0, ymax)
    plt.xlim(0.0, 1.0)
    
    plt.ylabel("Probability density",fontsize=fs)
    plt.xlabel(metric, fontsize=fs)
    plt.tick_params(axis='both', which='major', labelsize=fs-3)

    
    plt.scatter([0],[-1], color = 'm', alpha = 0.8, s=10,
                label= 'microbial (n='+str(len(microlist))+')')
        
    plt.scatter([0],[-1], color = 'c',alpha=0.8, s=10,
                label= 'macrobial (n='+str(len(macrolist))+')')
        
    plt.legend(bbox_to_anchor=(-0.03, 1.1, 3.75, .2), loc=10, ncol=2,
                        mode="expand",prop={'size':fs-2})
    
    
    
    
    
    fig.add_subplot(3, 3, 2)
    # scatter of Heips
    
    metric = 'Evenness'
    metlist = ObsHeip
        
    nmacs = 0
    nmics = 0
        
    indexList = range(len(klistSample))    
    #random.shuffle(indexList)
    
    for j in indexList:
        
        if klistSample[j] == 'm': 
            nmics += 1
            plt.scatter(NlistSample[j], metlist[j], color=klistSample[j], alpha=0.6, s=0.5)
        
        elif klistSample[j] == 'c':
            nmacs += 1
            plt.scatter(NlistSample[j], metlist[j], color=klistSample[j], alpha=0.6, s=0.5)
                                                                                                                                       
    plt.ylim(0.0, 1)
    plt.xscale('log')
    plt.xlim(1, 2*10**6)

    plt.xlabel('Total abundance', fontsize=fs)
    plt.ylabel(metric, fontsize=fs)
    plt.tick_params(axis='both', which='major', labelsize=fs-3)
    
    
    
    
    fig.add_subplot(3, 3, 3)
    """ Heips scatter feasible set """
    
    metric = 'Evenness'
    metlist = FS_Heip
        
    nmacs = 0
    nmics = 0
        
    for j in indexList:
        
        if klistSample[j] == 'm':
            nmics += 1
            plt.scatter(NlistSample[j], metlist[j], color=klistSample[j], alpha=0.6, s=0.5)
        
        elif klistSample[j] == 'c':
            nmacs += 1
            plt.scatter(NlistSample[j], metlist[j], color=klistSample[j], alpha=0.6, s=0.5)
        
                                                                                                                                            
    plt.ylim(0.6, 1)
    plt.xscale('log')
    plt.xlim(1, 2*10**6)

    plt.xlabel('Total abundance', fontsize=fs)
    plt.ylabel(metric, fontsize=fs)
    plt.tick_params(axis='both', which='major', labelsize=fs-3)
    
    
    
    
    fig.add_subplot(3, 3, 4)
    """ kdens of Berger-Parker """
    
    metric = 'Dominance'
    microlist = []
    macrolist = []
        
    for i, j in enumerate(BPlist):
        if klist[i] == 'm':
            microlist.append(j)
        elif klist[i] == 'c':
            macrolist.append(j)
            
    print len(microlist), len(macrolist)
                    
    D = get_kdens(microlist)
    n = len(microlist)
    plt.plot(D[0], D[1], color = 'm', lw=3, alpha = 0.9,
             label= 'microbial (n='+str(n)+')')
    ymax1 = max(D[1])

    
    n = len(macrolist)
    D = get_kdens(macrolist)
    plt.plot(D[0], D[1], color = 'c', lw=3, alpha = 0.9,
             label = 'macrobial (n='+str(n)+')')
    ymax2 = max(D[1])

    ymax = max(ymax1, ymax2) + 0.4*(max(ymax1, ymax2))
    plt.ylim(0.0, ymax)
    plt.ylabel("Probability density",fontsize=fs)
    plt.xlabel(metric, fontsize=fs)
    plt.tick_params(axis='both', which='major', labelsize=fs-3)    
    
    
    
    fig.add_subplot(3, 3, 5)
    # scatter of Berger-Parker
    
    metric = 'Dominance'
    metlist = ObsBP
        
    nmacs = 0
    nmics = 0
        
    indexList = range(len(klistSample))    
    #random.shuffle(indexList)
    
    for j in indexList:
        
        if klistSample[j] == 'm': 
            nmics += 1
            plt.scatter(NlistSample[j], metlist[j], color=klistSample[j], alpha=0.6, s=0.5)
            
        elif klistSample[j] == 'c':
            nmacs += 1
            plt.scatter(NlistSample[j], metlist[j], color=klistSample[j], alpha=0.6, s=0.5)
                
        
    plt.ylim(0, 1)
    plt.xscale('log')
    plt.xlim(10, 2*10**6)    
    
    plt.xlabel('Total abundance', fontsize=fs)
    plt.ylabel(metric, fontsize=fs)
    plt.tick_params(axis='both', which='major', labelsize=fs-3)
    
    
    
    fig.add_subplot(3, 3, 6)
    """ Berger-Parker scatter feasible set """
    
    metric = 'Dominance'
    metlist = FS_BP
        
    nmacs = 0
    nmics = 0
        
    for j in indexList:
        
        if klistSample[j] == 'm': 
            nmics += 1
            plt.scatter(NlistSample[j], metlist[j], color=klistSample[j], alpha=0.6, s=0.5)
        
        elif klistSample[j] == 'c':
            nmacs += 1
            plt.scatter(NlistSample[j], metlist[j], color=klistSample[j], alpha=0.6, s=0.5)
        
                                                                                                                                            
    plt.ylim(0, 0.5)
    plt.xscale('log')
    plt.xlim(10, 2*10**6)    

    plt.xlabel('Total abundance', fontsize=fs)
    plt.ylabel(metric, fontsize=fs)
    plt.tick_params(axis='both', which='major', labelsize=fs-3)
    
    
    
    
    fig.add_subplot(3, 3, 7)
    # kdens of % abundance of singletons
    
    metric = 'Number of singletons'
    microlist = []
    macrolist = []
        
    for i, j in enumerate(rareOnesList):
        if klist[i] == 'm':
            microlist.append(j)
        elif klist[i] == 'c':
            macrolist.append(j)
            
    print len(microlist), len(macrolist)
                
    D = get_kdens(microlist)
    n = len(microlist)
    plt.plot(D[0], D[1], color = 'm', lw=3, alpha = 0.9,
             label= 'microbial (n='+str(n)+')')
    ymax1 = max(D[1])

    
    n = len(macrolist)
    D = get_kdens(macrolist)
    plt.plot(D[0], D[1], color = 'c', lw=3, alpha = 0.9,
             label = 'macrobial (n='+str(n)+')')
    ymax2 = max(D[1])

    ymax = max(ymax1, ymax2) + 0.4*(max(ymax1, ymax2))
    
    plt.ylim(0, ymax)
    plt.xlim(1, 400)
    
    plt.ylabel("Probability density",fontsize=fs)
    plt.xlabel(metric, fontsize=fs)
    plt.tick_params(axis='both', which='major', labelsize=fs-3)
    
    
    
    fig.add_subplot(3, 3, 8)
    # scatter of % abundance of singletons
    
    metric = 'Number of singletons'
    metlist = ObsRare
        
    nmacs = 0
    nmics = 0
    
    indexList = range(len(klistSample))    
    #random.shuffle(indexList)
    
    for j in indexList:
        
        if klistSample[j] == 'm': 
            nmics += 1
            plt.scatter(NlistSample[j], metlist[j], color=klistSample[j], alpha=0.6, s=0.5)
        
        elif klistSample[j] == 'c':
            nmacs += 1
            plt.scatter(NlistSample[j], metlist[j], color=klistSample[j], alpha=0.6, s=0.5)
        
                
    plt.ylim(1, 10**3)
    plt.xlim(8, 2*10**6)
    plt.yscale('log')
    plt.xscale('log')

    plt.xlabel('Total abundance', fontsize=fs)
    plt.ylabel(metric, fontsize=fs)
    plt.tick_params(axis='both', which='major', labelsize=fs-3)
    
    
    
    fig.add_subplot(3, 3, 9)
    # scatter of % abundance of singletons
    
    metric = 'Number of singletons'
    metlist = FS_Rare
        
    nmacs = 0
    nmics = 0
    
    
    for j in indexList:
        
        if klistSample[j] == 'm': 
            nmics += 1
            plt.scatter(NlistSample[j], metlist[j], color=klistSample[j], alpha=0.6, s=0.5)
        
        elif klistSample[j] == 'c':
            nmacs += 1
            plt.scatter(NlistSample[j], metlist[j], color=klistSample[j], alpha=0.6, s=0.5)
        
                
    plt.ylim(1, 10**3)
    plt.xlim(8, 2*10**4)
    plt.yscale('log')
    plt.xscale('log')


    plt.xlabel('Total abundance', fontsize=fs)
    plt.ylabel(metric, fontsize=fs)
    plt.tick_params(axis='both', which='major', labelsize=fs-3)
    
    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/RareBio-Fig1.png', dpi=600)
    plt.close()
    
    
    return
    




#scatterMetrics()
#scatterMetricsRare()
#kdensMetrics()
#scatterNmax()
#scatterNmax2()

tool = 'Parts'
Fig1(tool)