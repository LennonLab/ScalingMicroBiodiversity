from __future__ import division
import sys
import os
import random
import numpy as np
from scipy import stats

mydir = os.path.expanduser("~/Desktop/Repos/rare-bio/")
sys.path.append(mydir + "tools/feasible_functions")
import feasible_functions as ff

mydir2 = os.path.expanduser("~/Desktop/")
OUT = open(mydir + 'output/EMPclosed-RADdata.txt','w+')

ct = 0
numEMP = 0
numEMPopen = 0
RADs = []
    
RADs = ff.GetSADsFromBiom_labeled(mydir2 +'data/micro', 'EMPclosed')

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
    SimpDom = ff.simpsons_dom(RAD)
            
    rareRel = ff.rarityRel(RAD)
    rareOnes = ff.rarityOnes(RAD)
            
    skew = float(np.log(stats.skew(RAD)))
    
    if S == 1: Var = 0
    elif S > 1: Var = np.var(np.array(RAD), ddof=1)
    
    else:
        print 'problem with Var'
    
    ct += 1
    
    print>>OUT, N, S, Evar, ESimp, ENee, EHeip, EQ, EPielou, BP, SimpDom, rareRel, rareOnes, skew
    print N, S, skew, len(RADs) - ct
    
print 'micros:', len(RADs)
    
OUT.close()
            
            
            
            
