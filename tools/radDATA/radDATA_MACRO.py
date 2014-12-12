from __future__ import division
import sys
import os
import random
import numpy as np
from scipy import stats

mydir = os.path.expanduser("~/Desktop/Repos/rare-bio")
sys.path.append(mydir + "/tools/feasible_functions")
import feasible_functions as ff

OUT = open(mydir + '/output/Macro-RADdata.txt','w+')
mydir2 = os.path.expanduser("~/Desktop")

ct = 0
numEMP = 0
numEMPopen = 0
RADs = []
    
    
def get_Macro_SADs():
    
    SADdict = {}
    for name in os.listdir(mydir2 +'/data/macro'):
        
        if name == 'BCI': continue
        if name == '.DS_Store': continue
            
        else:    
            print name      
            DATA = mydir2 + '/data/macro/'+name+'/'+name+'-data.txt'
                
            with open(DATA) as f: 
        
                for d in f:
                    if d.strip():
                        d = d.split()
                        
                        if name == 'GENTRY':
                            site = d[0]
                            #species = d[1] # Dataset name plus species identifier
                            abundance = float(d[2])
                        
                        else:
                            site = d[0]
                            #year = d[1]
                            #species = d[2] # Dataset name plus species identifier
                            abundance = float(d[3])
                
                        if abundance > 0:
                            if site in SADdict:
                                SADdict[site].append(abundance)
                	    else:
                                SADdict[site] = [abundance]
             
        
        
    RADs = []
    SADlist = SADdict.items()
    
    for tup in SADlist:
            
        RAD = tup[1]
        if len(RAD) >= 1: 
            RAD.sort()
            RAD.reverse()
            RADs.append(RAD)
            
    return RADs        
    
    
RADs = get_Macro_SADs()

for RAD in RADs:
        
    N = int(sum(RAD))
    S = len(RAD)
    
    if S < 1: continue 
    
    Evar = ff.e_var(RAD)
    ESimp = ff.simpsons_evenness(RAD)
    ENee = ff.NHC_evenness(RAD)
    EPielou = ff.pielous_evenness(RAD)
            
    EHeip = ff.Heips_evenness(RAD)
    EQ = ff.EQ_evenness(RAD)
            
    BP = ff.Berger_Parker(RAD)
    SimpDom = 0 #ff.simpsons_dom(RAD)
            
    rareRel = ff.rarityRel(RAD)
    rareOnes = ff.rarityOnes(RAD)
            
    skew = float(np.log(stats.skew(RAD)))
    
    if S == 1: Var = 0
    elif S > 1: Var = np.var(np.array(RAD), ddof=1)
    
    else:
        print 'problem with Var'
    
    ct += 1
    
    print>>OUT, N, S, Evar, ESimp, ENee, EHeip, EQ, EPielou, BP, SimpDom, rareRel, rareOnes, skew
    #            N, S, ESimp, EHeip, BP, SimpDom, rareRel, rareOnes, rareSumOnes, Var = data_list
        
    print N, S, skew, len(RADs) - ct
    
print 'macros:', len(RADs)
    
OUT.close()
            
            
            
            
