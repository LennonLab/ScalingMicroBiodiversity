import sys

DATA = '/Users/lisalocey/Desktop/RareBio/data/emp/emp/emp-SSADdata.txt'
    
SADdict = {}
    
with open(DATA) as f: 
    
    for d in f:
        if d.strip():
            
            print d
            d = d.split()
            print d
            
            sample = d[1]
            abundance = float(d[2])
            
            print sample, abundance
            sys.exit()
            if abundance > 0:
                if sample not in SADdict:
                    SADdict[sample] = [abundance]
        
                else: SADdict[sample].append(abundance)
    