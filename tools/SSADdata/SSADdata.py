from __future__ import division
import sys
import os
import random
import numpy as np

sys.path.append("/Users/lisalocey/Desktop/Repos/rare-bio/tools/feasible_functions")
import feasible_functions as ff

mydir = os.path.expanduser("~/Desktop")
OUT = open(mydir+'/Repos/rare-bio/output/EMPclosed-SSAD-ResultsTable.txt','w+')


def get_SSADs():

    DATA = mydir + '/data/micro/EMP/EMPclosed/EMPclosed-SSADdata.txt'
    
    SADdict = {}
    
    with open(DATA) as f: 
        
        for d in f:
            if d.strip():
                
                d = d.split()
                species = d[0]
                #sample = d[1]
                abundance = float(d[2])
                
                if abundance > 0:
                    if species in SADdict:
                        SADdict[species].append(abundance)
                    else:
                        SADdict[species] = [abundance]
             
        
        
    SADs = []
    SADlist = SADdict.items()
    
    for tup in SADlist:
            
        SAD = tup[1]
        if len(SAD) >= 1: 
            SAD.sort()
            SAD.reverse()
            SADs.append(SAD)
            
    return SADs        
    
    



ct = 0
numEMP = 0
numEMPopen = 0

SSADs = get_SSADs()
num = len(SSADs)

for SSAD in SSADs:
    
    #SSAD = list([x for x in SSAD if x != 0]) # removes zeros
    
    N = int(sum(SSAD))
    S = int(len(SSAD))
    print N, S,'  ',num - ct
            
    if S < 10: continue 
    
    Evar = ff.e_var(SSAD)
    ESimp = ff.simpsons_evenness(SSAD)
    ENee = ff.NHC_evenness(SSAD)
    EPielou = ff.pielous_evenness(SSAD)
    
    EHeip = ff.Heips_evenness(SSAD)
    EQ = ff.EQ_evenness(SSAD)
    
    BP = ff.Berger_Parker(SSAD)
    SimpDom = ff.simpsons_dom(SSAD)
    perc_ones = ff.Singletons(SSAD)
    
    SSAD.sort()
    SSAD.reverse()
    
    rareRel = ff.rarityRel(SSAD)
    rareOnes = ff.rarityOnes(SSAD)
    rarePair = ff.rarityPair(SSAD)
    rareSumOnes = ff.raritySumOnes(SSAD)
    ct+=1
    
    print>>OUT, N, S, Evar, ESimp, ENee, EHeip, EQ, EPielou, BP, SimpDom, rareRel, rareOnes, rarePair, rareSumOnes 
    
print 'number of SSADs:', len(SSADs)

OUT.close()