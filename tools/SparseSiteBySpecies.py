import os

mydir2 = os.path.expanduser("~/Desktop/")

OrC = 'open'

def GetSADsFromBiom_labeled(path, dataset):
    
    minS = 2
    
    IN = path + '/' + dataset + '-SSADdata.txt'
    n = sum(1 for line in open(IN))
    
    SADdict = {}
    
    with open(IN) as f: 
        
        for d in f:
            
            print 'Reading in SSAD data. Lines left:', n
            n -= 1
            
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
    
    IN.close()
    
    OUT = open(path + '/' + dataset + '-SbyS.txt','w+')
        
    SADlist = SADdict.items()
    n = len(SADlist)
    
    for i, tup in enumerate(SADlist):
                
        SAD = tup[1]
        if len(SAD) >= minS: 
            SAD.sort()
            SAD.reverse()
            print n - i
            print >> OUT, SAD
                
                    
    OUT.close()
        
    return 
        
GetSADsFromBiom_labeled(mydir2 +'data/micro/EMP'+OrC, 'EMP'+OrC)