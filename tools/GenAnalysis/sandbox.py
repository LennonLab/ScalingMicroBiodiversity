import random

#print random.randint(0,1)
#print range(1,10)

mlist = [2]

#print random.randint(0,1)


N = 5*(10**30)
P = 2.9*10**27
S = 7.0*10**26

#print P/N*100, S/N*100

def RandomComposition(Q,N):
    composition = []
    
    indices = random.sample(range(1,Q+1),N-1) # get a random sample of k-1 objects from n-1 objects
    indices.sort()
    Qlist = [1]*Q
    
    partsOF = [Qlist[i:j] for i, j in zip([0]+indices, indices+[None])]
    for p in partsOF:
        composition.append(sum(p))
    
    return composition
    

def RanComp(Q, N):
    composition = []
    
    indices = []
    while len(indices) < N-1:
        index = random.randint(1, Q)
        if index in indices: continue
        else: indices.append(index)
    
    indices.sort()
    Qlist = [1]*Q
    
    partsOF = [Qlist[i:j] for i, j in zip([0]+indices, indices+[None])]
    for p in partsOF:
        composition.append(sum(p))
    
    return composition
        
                
RAD = RandomComposition(10**4, 10**1)
print max(RAD)