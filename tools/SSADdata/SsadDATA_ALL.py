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



ct = 0
numMicros = 0
numMacros = 0
numEMP = 0
numEMPopen = 0
    
for dataset in datasets:
    
    name = dataset[0] # name of dataset
    if name == '.DS_Store': continue
    
    kind = dataset[1] # micro or macro
    OUT = open(mydir2 + 'data/'+kind+'/'+name+'/'+name+'-SSADMetricData.txt','w+')
    RADs = []
    
    if kind == 'macro':
        path = kind
        RADs = ff.get_SADs(mydir2 +'/data/'+path, name)

        print 'macro', name, len(RADs)
        
        
    if kind == 'micro':
        #RADs = getPhyloSADs.get_SADs(seqID, dataset[0], 'genus')    

        path = kind
        if name == 'EMPclosed': continue
        if name == 'EMPopen':
            RADs = ff.GetSADsFromBiom_labeled(mydir2 +'data/'+path+'/EMPopen/', name)
        
        else: RADs = ff.get_SADs(mydir2 +'data/'+path, name)
                                                                                                        
    ct = 0
    numRADs = len(RADs)
    for RAD in RADs:
        
        RAD = list([x for x in RAD if x != 0])
        
        N = sum(RAD)
        S = len(RAD)
        
        if S < 2: continue 
        
        Evar = ff.e_var(RAD)
        ESimp = ff.simpsons_evenness(RAD)
        ENee = ff.NHC_evenness(RAD)
        EPielou = ff.pielous_evenness(RAD)
        
        EHeip = ff.Heips_evenness(RAD)
        EQ = ff.EQ_evenness(RAD)
        
        BP = ff.Berger_Parker(RAD)
        SimpDom = ff.simpsons_dom(RAD)
        perc_ones = ff.Singletons(RAD)
        
        Nmax = max(RAD)
    
        rareRel = ff.rarityRel(RAD)
        rareOnes = ff.rarityOnes(RAD)
        
        skew = stats.skew(RAD)
        
        ct+=1
        
        print>>OUT, name, kind, N, S, Evar, ESimp, ENee, EHeip, EQ, EPielou, BP, SimpDom, rareRel, rareOnes, skew
        print 'micro', name, numRADs - ct
        
    OUT.close()