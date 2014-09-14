from __future__ import division
import sys

import matplotlib.pyplot as plt
import feasible_functions as ff
import random
from random import randrange, randint, uniform, seed, choice

import time
import numpy as np
import scipy.stats                                       

    


def RandCompFast(Q, N):
    
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

    return composition





def SimLogNorm(N,S, sample_size=1):
    '''This script codes the Lognormal Model'''
    sample = []
    while len(sample) < sample_size:
        
        n = 0.75 * N
        RAC = [n, N - n]
        
        while len(RAC) < S:
            
            ind = randrange(len(RAC))
            v = RAC.pop(ind)
            v1 = round(0.75 * v) 
            v2 = v - v1   # forcing all abundance values to be integers
            
            RAC.extend([v1, v2])
        
        if len(RAC) > 1:
            RAC.sort()
            RAC.reverse()
            sample.append(RAC)
    
    return sample




def DomPreInt(N, sample_size, rel=False): # Works only with positive integers
    '''This script codes Tokeshi's Dominance Preemption Model
    this code does not work well with small N or high S'''
    sample = [] # A list of RACs
   
    while len(sample) < sample_size: # number of samples loop     
        RAC = [] #RAC is a list
        sp1 = randrange(int(round(N *.5)), N) #Rand select from N to N(.5)
        ab1 = N - sp1
        RAC.extend([sp1, ab1])
        
        while RAC:
            
            ab2 = RAC.pop()
            if ab2 <= 1:
                break
                
            sp2 = randrange(int(round(ab2*.5)), ab2)
            RAC.extend([sp2, ab2-sp2])
        
        if min(RAC) < 1:
            print 'min',min(RAC)
            sys.exit()
            
        if len(RAC) > 2: sample.append(RAC)
    return sample





def SimParetoInt(N, sample_size, rel=False):
    '''This script codes the Pareto Model'''
    sample = []
    while len(sample) < sample_size: 
        RAC = [0.8*N, 0.2*N]
        
        while RAC:
            ind = randrange(len(RAC))
            v = RAC.pop(ind)
            v1 = int(round(0.8 * v))
            v2 = v - v1  # forcing all abundance values to be integers
            
            if v1 < 1 or v2 < 1: break  # forcing smallest abundance to be 
                                        # greater than one
            RAC.extend([v1, v2])       
                        
        if len(RAC) > 1:        
            RAC.sort(reverse = True)
            sample.append(RAC)
    
    return sample
    
    
def Sample_SimpleRandomFraction(N, S, sample_size):
    
    """ 
    This function randomly and sequently splits N into two integers by
    starting with a vector where N is the only value, i.e. N, and then splitting
    N into two positive integers at random, [N] -> [x1, x2], where x1+x2 = N.
    The function then draws one of the integers at random from the new vector
    (having 2 values), and then splits the integer into two other integers. 
    At this point, the vector has three values [x1, x2, x3], where x1+x2+x3 = N.  
    This process keeps on until there are a number of integers equal to S.
    
    N  :  total abundance; number of individuals in the community
    S  :  species richness; number of species in the community
    sample_size  :  number of random rank-abundance vectors to generate    
    """
    
    sample = []
    for i in range(sample_size):
        
        RAC = [N]
        while len(RAC) < S:
                
            sp1_ab = choice(RAC) # pick a species (i.e. vector value) abundance at random
            if sp1_ab == 1:
                continue # this is a control statement to prevent the program
                # from encountering an impossible conditions, i.e., dividing 1
                # at random into two positive integers
                
            sp2_ab = randrange(1, sp1_ab) # pick a random number (abundance) between 1 and sp1_ab - 1, inclusive
            sp1_index = RAC.index(sp1_ab) # index in the RAC of the species we picked
                
            RAC[sp1_index] = sp1_ab - sp2_ab # decrease its abundance according to sp_ab
            RAC.append(sp2_ab)
            
        RAC.sort(reverse = True)  #sorts RAC's in decending order to aid in graphing. 
        sample.append(RAC) # appending a list (i.e. an RAC) to another list
    
    return sample
    
    



def DomDecayInt(N, sample_size, rel=False): # Works only with positive integers
    sample = [] # A list of RACs
   
    while len(sample) < sample_size: # number of samples loop     
        
        RAC = [] #RAC is a list
        sp1 = randint(1, int(round(N*.5)))
        ab1 = N - sp1
        RAC.extend([sp1, ab1]) 
        
        while RAC:
            if min(RAC) < 2: break
            
            ab2 = RAC.pop()
            sp2 = randint(1, int(round(ab2 * .5)))
            RAC.extend([sp2, ab2 - sp2])

        if len(RAC) > 1:
            RAC.sort(reverse = True)
            sample.append(RAC)
     
    return sample




def RADplot():
    
    N = 10**6 #sum(GSAD)
    S = 10**3 #len(GSAD)
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    plt.title('Total abundance = '+str('%.1e' % N)+', Richness = '+str('%.1e' % S))
    
    """
    # Observed
    print 'Berger Parker as % =', 100*max(GSAD)/N,'\n'
    ranks = range(1, S+1)
    plt.scatter(ranks, GSAD, lw = 0, c='k', label='Pooled EMP', alpha=0.5) 
    """
    
    # Predicted geometric series, Ken
    t0 = time.time()
    predRAD = RandCompFast(N, S)
    t = time.time() - t0
    print 'time for Geometric series:',t
    print 'Berger Parker as % =', 100*max(predRAD)/sum(predRAD),'\n'
    ranks = range(1, len(predRAD)+1)
    plt.plot(ranks, predRAD, lw = 1, c='m', label='Geometric Series') 
    
    """
    # Dominance Preemption
    t0 = time.time()
    DPlist = DomPreInt(N, 1)
    predRAD = DPlist[0]
    t = time.time() - t0
    print 'time for Dominance Preemption:',t
    print 'Berger Parker as % =', 100*max(predRAD)/sum(predRAD),'\n'
    ranks = range(1, len(predRAD)+1)
    plt.plot(ranks, predRAD, lw = 1, c='c', label='Dominance Preempt.') 
    
    # Dominance Decay
    t0 = time.time()
    DDlist = DomDecayInt(N, 1)
    predRAD = DDlist[0]
    t = time.time() - t0
    print 'time for Dominance Decay:',t
    print 'Berger Parker as % =', 100*max(predRAD)/sum(predRAD),'\n'
    ranks = range(1, len(predRAD)+1)
    plt.plot(ranks, predRAD, lw = 1, c='Steelblue', label='Dominance Decay') 
    """
    
    
    # Sim Log-normal
    t0 = time.time()
    LNlist = SimLogNorm(N, S)
    predRAD = LNlist[0]
    t = time.time() - t0
    print 'time for Log-normal:',t
    print 'Berger Parker as % =', 100*max(predRAD)/sum(predRAD),'\n'
    ranks = range(1, len(predRAD)+1)
    plt.plot(ranks, predRAD, lw = 1, c='r', label='Log-normal', alpha=0.8) 
    
    
    # Simple RandFrac
    t0 = time.time()
    LNlist = Sample_SimpleRandomFraction(N, S, 1)
    predRAD = LNlist[0]
    t = time.time() - t0
    print 'time for RandFrac:',t
    print 'Berger Parker as % =', 100*max(predRAD)/sum(predRAD),'\n'
    ranks = range(1, len(predRAD)+1)
    plt.plot(ranks, predRAD, lw = 1, c='c', label='RandFrac', alpha=0.8) 
    
    
    leg = plt.legend(loc=3,prop={'size':16})
    leg.draw_frame(False)
    
    plt.xlabel('Rank in abundance', fontsize = 16)
    plt.ylabel('Abundance', fontsize = 16)
    plt.yscale('log')
    #plt.xscale('log')
    
    plt.savefig('/Users/lisalocey/Desktop/RareBio/global/figs/All_EMP_N='+str(int(N))+'_S='+str(int(S))+'.png',dpi=600)
    plt.show()
    
    return







SADs = []    
SpDict = {}

"""
SADdict = ff.GetSADsFromBiom_labeled('/Users/lisalocey/Desktop/RareBio/data/micro/EMP', 'EMP')
SADlist = SADdict.items()
            
for SAD in SADlist:
    site = SAD[0]
    lRAD = SAD[1]
    
    
    for _list in lRAD:
        sp = _list[0]
        ab = _list[1]
        
        if sp in SpDict:
            SpDict[sp] = SpDict[sp] + ab
        
        else: SpDict[sp] = ab
            

GSAD = []            
DSAD = SpDict.items()

for _list in DSAD:
    GSAD.append(_list[1])
    

print 'N = ', sum(GSAD)
print 'S = ', len(GSAD)

#GSAD = random.sample(GSAD, 100)
#GSAD = np.random.logseries(0.97, 10^3)  # list of propagules
#GSAD = GSAD.tolist()
GSAD.sort()
GSAD.reverse()
"""

RADplot()
    


