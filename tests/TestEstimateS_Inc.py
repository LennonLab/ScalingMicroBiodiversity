from __future__ import division
import  matplotlib.pyplot as plt
import numpy as np
import random
import os
import sys

#import scipy
from scipy import stats
from random import randint

mydir = os.path.expanduser("~/Desktop/Repos/rare-bio/tools/StatPak")
sys.path.append(mydir)
import DiversityMetrics as DM

def Chao2(SbySdict):

    SpList = []
    for key in SbySdict:
        SpList.extend(SbySdict[key])
        
    SpList = list(set(SpList))
    S = len(SpList)
    
    SpIncList = [0]*len(SpList)
    Sites = [0]*len(SbySdict.keys())   
    
    for i, sp in enumerate(SpList):
        print i, sp, len(SpList)-i
        for key in SbySdict:
            if sp in SbySdict[key]:
                SpIncList[i] += 1
    
    for i in SpIncList:
        Sites[i-1] += 1
    
    q1 = Sites[0]
    q2 = Sites[1]
    m = sum(Sites)

    chao2 = S + (((m-1)/m) * ((q1*(q1-1)) / (2*(q2+1))))

    return [int(chao2), S]


# Generate Site by Species matrix
SbyS = []
IN = '/Users/lisalocey/Desktop/data/micro/EMPclosed/EMPclosed-SSADdata.txt'
SbySdict = {}
    
with open(IN) as f:
        
    for d in f:
        if d.strip():
            
            d = d.split()
            species = d[0]
            sample = d[1]
            abundance = float(d[2])
            
            if abundance > 0:
                if sample not in SbySdict:
                    SbySdict[sample] = [species]
                    
                else: SbySdict[sample].append(species)


obsSample = SbySdict.keys()
n = len(obsSample)


pdist = []
print 'Number of sites:', n

Chao2s = []
uChao2SEs = []
lChao2SEs = []

Ns = []
Ss = []

n = 10
while n < 1000:
    
    for i in range(4):
        sC = []
        ns = []
        ss = []
        
        SiteList = np.random.choice(obsSample, size=n, replace=True)
        SbS_trunc = {}
        
        for name in SiteList: 
            SbS_trunc[name] = list(SbySdict[name])
        
        chao2, S = Chao2(SbS_trunc)
        sC.append(chao2)
        
    print 'Sample Size:',n, ' Richness:',S,' Chao2:', chao2
    #sys.exit()
        
    avgCh2 = float(np.mean(sC))
    Chao2s.append(avgCh2)
    lChao2SEs.append(avgCh2 - float(stats.sem(sC)))
    uChao2SEs.append(avgCh2 + float(stats.sem(sC)))
    
    print n
    Ns.append(n)
    n = n*2
    

fig = plt.figure()
fig.add_subplot(1,1,1)

plt.plot(Ns, Chao2s, color='r', alpha=0.8)
#plt.plot(Ns, uChao1SEs, color='r', ls='--', alpha=0.6)
#plt.plot(Ns, lChao1SEs, color='r', ls='--', alpha=0.6)

plt.plot(Ns, Ss, color='0.3', ls='-')

plt.ylabel('Richness')
plt.xlabel('Samples')

plt.xscale('log')
plt.show()