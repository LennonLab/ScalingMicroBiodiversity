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
import Chao1, ACE


def Chao2(SbySdict):

    obsSpecies = SbySdict.keys()
    S = len(obsSpecies)
    fs = [0]*S    
    
    counts = []
    Zs = False

    for sp in obsSpecies:
        f = len(SbySdict[sp])
        counts.append(f)
        fs[f] += 1
    
    if Zs is False: fs.pop(0)
    
    q1 = fs[0]
    q2 = fs[1]
    m = sum(fs)

    chao2 = S + (((m-1)/m) * ((q1*(q1-1)) / (2*(q2+1))))

    return [int(chao2), counts]


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
                if species not in SbySdict:
                    SbySdict[species] = [abundance]
                    
                else: SbySdict[species].append(abundance)


obsSpecies = SbySdict.keys()
S = len(obsSpecies)

gSAD = []
for sp in obsSpecies:
    gSAD.append(int(sum(SbySdict[sp])))

N = sum(gSAD)

pdist = []
for i, val in enumerate(gSAD): pdist.append(val/N)

print 'Global N:', N
print 'Global S:', S
print 'sum pdist:',sum(pdist)

Chao1s = []
uChao1SEs = []
lChao1SEs = []

ACEs = []
uACESEs = []
lACESEs = []

Ns = []
Ss = []

n = 1000
while n < 1000:
    for i in range(100):
        sC = []
        sA = []
        ns = []
        ss = []
        
        if len(obsSpecies) != len(pdist):
            print 'len(obsSpecies) != len(pdist)'
            sys.exit()
        
        abv = np.random.choice(obsSpecies, size=n, replace=True, p=pdist)
        print abv
        abv = abv.tolist()
        print abv
        
        splist = list(set(abv))
        print splist, len(splist)
        
        ss.append(len(splist))
        
        sad = []
        for sp in splist: sad.append(abv.count(sp))
        
        sad.sort()
        sad.reverse()
        print sad
        sys.exit()
        
        
        ch1 = Chao1.chao1(sad, bias_corrected=True)
        
        print ch1, len(sad), sum(sad), n, len(abv)
        
        
        sC.append(ch1)
        Ace = ACE.ace(sad, rare_threshold=10)
        sA.append(Ace)
    
    avgCh1 = float(np.mean(sC))
    Chao1s.append(avgCh1)
    lChao1SEs.append(avgCh1 - float(stats.sem(sC)))
    uChao1SEs.append(avgCh1 + float(stats.sem(sC)))
    
    avgAce = float(np.mean(sA))
    ACEs.append(avgAce)
    lACESEs.append(avgAce - float(stats.sem(sA)))
    uACESEs.append(avgAce + float(stats.sem(sA)))
    
    print n
    Ss.append(np.mean(ss))
    Ns.append(n)
    
    if np.mean(ss) > n:
        print 'S > N', np.mean(ss), n
        sys.exit()
        
    n = n*2
    

fig = plt.figure()
fig.add_subplot(1,1,1)

plt.plot(Ns, Chao1s, color='r', alpha=0.8)
#plt.plot(Ns, uChao1SEs, color='r', ls='--', alpha=0.6)
#plt.plot(Ns, lChao1SEs, color='r', ls='--', alpha=0.6)

plt.plot(Ns, ACEs, color='b', alpha=0.8)
#plt.plot(Ns, uACESEs, color='b', ls='--', alpha=0.6)
#plt.plot(Ns, lACESEs, color='b', ls='--', alpha=0.6)

plt.plot(Ns, Ss, color='0.3', ls='-')

plt.ylabel('Richness')
plt.xlabel('Individuals')

plt.xscale('log')
plt.show()
   

#print '\nBased on incidence:'
#chao2, counts = Chao2(SbyS)
#print 'Estimated via Chao2:', chao2#, 'Chao1 Confidence intervals:', chao1CI
