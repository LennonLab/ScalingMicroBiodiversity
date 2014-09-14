from __future__ import division
import sys                                            
import numpy as np
import math
import random
import matplotlib.pyplot as plt
import scipy.stats                                       
from random import randrange, randint, uniform, seed, choice

    

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



    
def SimLogNormInt(N, sample_size, rel=False):
    '''This script codes the Lognormal Model'''
    sample = []
    while len(sample) < sample_size:
        
        n = int(round(0.75 * N))
        RAC = [n, N - n]
        
        while RAC:
            
            ind = randrange(len(RAC))
            v = RAC.pop(ind)
            v1 = int(round(0.75 * v)) 
            v2 = v - v1   # forcing all abundance values to be integers
            
            if v1 < 1 or v2 < 1: break  # forcing smallest abundance to be 
                                        # greater than one
            RAC.extend([v1, v2])
        
        if len(RAC) > 1:
            RAC.sort()
            RAC.reverse()
            sample.append(RAC)
    
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
    
    
def Sample_SimpleRandomFraction(N, S, sample_size, rel=False):
    
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
        
        while RAC:
                
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

