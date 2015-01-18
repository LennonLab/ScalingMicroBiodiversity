import sys
import os

mydir = os.path.expanduser("~/Desktop/Repos/rare-bio/")
mydir2 = os.path.expanduser("~/Desktop/")

OrC = 'open'

def GetSADsFromBiom_labeled(path, dataset):
    
    minS = 2
    
    IN = path + '/' + dataset + '-SSADdata.txt'
    SADdict = {}
    
    with open(IN) as f: 
        
        for d in f:
            if d.strip():
                
                d = d.split()
                species = d[0]
                sample = d[1]
                abundance = float(d[2])
                
                if abundance > 0:
                    if sample not in SADdict:
                        
                        SADdict[sample] = [[species, abundance]]
            
                    else:
                        SADdict[sample].append([species, abundance])
    
    OUT = open(mydir + 'output/EMP'+OrC+'-SbyS.txt','w+')
        
    SADlist = SADdict.items()
    
    for tup in SADlist:
                
        SAD = tup[1]
        if len(SAD) >= minS: 
            SAD.sort()
            SAD.reverse()
            print >> OUT, SAD
                
                    
    OUT.close()
        
    return 
        
GetSADsFromBiom_labeled(mydir2 +'data/micro/EMP'+OrC, 'EMP'+OrC)