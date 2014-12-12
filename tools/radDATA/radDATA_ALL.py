from __future__ import division
import sys
import os
import random
import numpy as np
from scipy import stats

mydir = os.path.expanduser("~/Desktop/Repos/rare-bio/")
mydir2 = os.path.expanduser("~/Desktop/")


sys.path.append(mydir + "tools/feasible_functions")
import feasible_functions as ff



def rareTail(RAD):
    
    RAD.sort() # RAD ordered small to large
    N = float(sum(RAD))
    S = float(len(RAD))
    
    eighty = 0.8*N
    ten = 0.1*N
    one = 0.01*N
    ptone = 0.001*N
    ptzone = 0.0001*N
    
    t = 0
    s1 = 0
    s2 = 0
    s3 = 0
    s4 = 0
    s5 = 0
    
    for ab in RAD:
        t += ab
        if t <= ptzone:
            s1 += 1
        if t <= ptone:
            s2 += 1
        if t <= one:
            s3 += 1
        if t <= ten:
            s4 += 1
        if t <= eighty:
            s5 += 1
        if t > eighty: break
        #t = 0
            
    return [s5/S, s4/S, s3/S, s2/S, s1/S]
    

def domTail(RAD):
    
    RAD.reverse()
    def getRel(num, RAD):
        ct = 0
        t = 0
        while ct < num:
            t += RAD[ct]
            ct+=1 
        return t
    
    rel = getRel(1, RAD)    
    BP1 = rel/sum(RAD)
    rel = getRel(2, RAD)
    BP2 = rel/sum(RAD)
    rel = getRel(3, RAD)
    BP3 = rel/sum(RAD)   
    rel = getRel(4, RAD)
    BP4 = rel/sum(RAD)
    rel = getRel(5, RAD)
    BP5 = rel/sum(RAD)

    return [BP1, BP2, BP3, BP4, BP5]



datasets = []
seqID = '99'

"""
for name in os.listdir('/Users/lisalocey/Desktop/RareBio/data/micro9599/'+seqID+'data'): 
    index = name.find('-')
    name = name[0:index]
    datasets.append([name, 'micro'])        
"""

for name in os.listdir(mydir2 +'data/micro'):
        datasets.append([name, 'micro'])   

for name in os.listdir(mydir2 +'data/macro'):
        datasets.append([name, 'macro'])   


OUT = open(mydir + 'output/RADdata.txt','w+')
#OUT = open(mydir + 'output/Macro-RADdata.txt','w+')


ct = 0
numMicros = 0
numMacros = 0
numEMP = 0
numEMPopen = 0
    
for dataset in datasets:
    
    name = dataset[0] # name of dataset
    kind = dataset[1] # micro or macro
    RADs = []
    
    if name == '.DS_Store': continue
    
    
    if kind == 'micro':
        #RADs = getPhyloSADs.get_SADs(seqID, dataset[0], 'genus')    

        path = kind
        
        if name == 'EMPclosed': # Use open or closed OTU reference
            #continue
            RADs = ff.GetSADsFromBiom_labeled(mydir2 +'data/'+path+'/EMPopen/', name)
        
        if name == 'EMPopen':
            continue
            RADs = ff.GetSADsFromBiom_labeled(mydir2 +'data/'+path+'/EMPopen/', name)
        
        else: RADs = ff.get_SADs(mydir2 +'/data/'+path, name)
        
        print 'micro', len(RADs)
        
        if len(RADs) > 200:
            RADs = random.sample(RADs, 180) # getting between 180 and 200 SADs per dataset    
        
        numMicros += len(RADs)
                                        
                                                                                                        
    if kind == 'macro':
        path = kind
        RADs = ff.get_SADs(mydir2 +'/data/'+path, name)

        print 'macro', len(RADs)
        
        if len(RADs) > 200:
            RADs = random.sample(RADs, 100)
        
        numMacros += len(RADs)

    
    
    for RAD in RADs:
        
        RAD = list([x for x in RAD if x != 0])
        
        N = sum(RAD)
        S = len(RAD)
        
        if S < 10: continue 
        
        Evar = ff.e_var(RAD)
        ESimp = ff.simpsons_evenness(RAD)
        ENee = ff.NHC_evenness(RAD)
        EPielou = ff.pielous_evenness(RAD)
        
        EHeip = ff.Heips_evenness(RAD)
        EQ = ff.EQ_evenness(RAD)
        
        BP = ff.Berger_Parker(RAD)
        SimpDom = ff.simpsons_dom(RAD)
        perc_ones = ff.Singletons(RAD)
        
        Nmax = np.percentile(RAD, 100)
    
        rareRel = ff.rarityRel(RAD)
        rareOnes = ff.rarityOnes(RAD)
        
        skew = stats.skew(RAD)
        
        ct+=1
        
        print>>OUT, name, kind, N, S, Evar, ESimp, ENee, EHeip, EQ, EPielou, BP, SimpDom, rareRel, rareOnes, skew
        #print>>OUT,             N, S, Evar, ESimp, ENee, EHeip, EQ, EPielou, BP, SimpDom, rareRel, rareOnes, skew
    
print 'micros:',numMicros, ' macros:',numMacros
    
OUT.close()
            
            
            
            
