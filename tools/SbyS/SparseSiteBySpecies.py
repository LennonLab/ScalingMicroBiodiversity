import os

mydir2 = os.path.expanduser("~/Desktop/")

OrC = 'closed'

def GetSADsFromBiom_labeled(path, dataset):
    
    minS = 2
    
    IN = path + '/' + dataset + '-SSADdata.txt'
    n = sum(1 for line in open(IN))
    
    SiteDict = {}
    
    print 'Starting build of SiteDict'
    with open(IN) as f: 
        
        for d in f:
            
            #print 'Reading in SSAD data. Lines left:', n
            #n -= 1
            
            if d.strip():
                
                d = d.split()
                species = d[0]
                sample = d[1]
                abundance = float(d[2])
                
                if abundance > 0:
                    if sample not in SiteDict:
                        
                        SiteDict[sample] = [species]
            
                    else:
                        SiteDict[sample].append(species)
    
    print 'Finished building SiteDict'
    
    
    OUT = open(path + '/' + dataset + '-SbyS.txt','w+')
        
    SiteLists = SiteDict.items()
    n = len(SiteLists)
    
    for i, site in enumerate(SiteLists):
                
        if len(site) >= minS: 
            print n - i
            print >> OUT, site
                    
    OUT.close()
    
    print 'Finished generating SbyS file'    
    return 
        
GetSADsFromBiom_labeled(mydir2 +'data/micro/EMP'+OrC, 'EMP'+OrC)