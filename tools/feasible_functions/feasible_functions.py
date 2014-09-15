from __future__ import division
import sys
import numpy as np
from scipy import stats
from scipy.stats import gaussian_kde
import re
import random

"""
sys.path.append("/Users/lisalocey/Desktop/RareBio/")
import ModelsN

sys.path.append("/Users/lisalocey/Desktop/RareBio/global/GenAnalysis/tools/")
import mete
#import pln
"""

sys.path.append("/Users/lisalocey/Desktop/Repos/partitions/partitions/Python/")
import partitions


########################################################################################################
######   A Section devoted to evenness indices and descriptive statistical functions ###################
minS = 7


def rarityRel(sad): # relative abundance of the least abundant taxon
    
    """ the probability that the next individual sampled will be from the least 
    abundant taxon """
    
    return min(sad)/sum(sad)
    
    
def rarityOnes(sad): 
    
    """ The rarest you can be is to have a single member. The more singletons,
    the more rarity """
    
    return sad.count(1)

    
def rarityPair(sad): 
    """ fraction of times (probability) that the next sampled individual comes
    from a taxon of lesser abundance, in a pair of ordered samples """
    
    weights = np.array(sad)/sum(sad)
    n = 1000
    ct = 0
    f = 0
    
    while ct < n:
        s = np.random.choice(sad, size=1, replace=True, p=weights)
        s1 = np.random.choice(sad, size=1, replace=True, p=weights)
        if s1 < s:
            f += 1
        
        ct+=1
        
    return f/n
    

    

def raritySumOnes(sad):
    """ The rarest you can be is to have a single member. The more singletons,
    the more rarity. """
    
    return sad.count(1)/sum(sad)


    

def Berger_Parker(sad):
    return max(sad)/sum(sad)

def Singletons(sad):
    singletons = sad.count(1)
    return 100*(singletons/len(sad))

def Shannons_H(sad):
    H = 0
    for i in sad:
        p = i/sum(sad)
        H += p*np.log(p)
    return H*-1.0
    
def Shannons_even(sad):
    H = Shannons_H(sad)
    S = len(sad)
    return H/np.log(S)
    
def simplest_gini(x):
        """Return computed Gini coefficient of inequality. 
        This function was found at http://econpy.googlecode.com/svn/trunk/pytrix/utilities.py """

        #note: follows basic formula
        #see: `calc_gini2`
        #contact: aisaac AT american.edu
        	
        x = sorted(x)  # increasing order
        n = len(x)
        G = sum(xi * (i+1) for i,xi in enumerate(x))
        G = 2.0*G/(n*sum(x)) #2*B
        return G - 1 - (1./n)

def gini_sample(SADs):
    """ Compute Gini's coefficient for each macrostate in a random sample """
    Gs = []
    for sad in SADs:
        G = simplest_gini(sad)
        Gs.append(G)
    return Gs


def Mcintosh_evenness(SAD):
    S = len(SAD)
    N = sum(SAD)
    sum_n = 0
    for n in SAD: sum_n += n**2
    U = np.sqrt(sum_n)    
    E = (N - U)/(N - (N/np.sqrt(S)))
    return E


def pielous_evenness(SAD):
    S = len(SAD)
    N = float(sum(SAD))
    H = 0
    for p in SAD:
        H += -(p/N)*np.log(p/N)
    J = H/np.log(S)
    return J


def NHC_evenness(SAD):
    SAD.sort()
    SAD.reverse()
    x_list = range(1,len(SAD)+1)
    y_list = np.log(SAD)
    slope,intercept,r_value,p_value,std_err = stats.linregress(x_list, y_list)
    
    if slope > 0.0:
        evar = e_var(SAD)
        print slope, p_value, evar
    return slope


def Heips_evenness(SAD):
    S = len(SAD)
    N = float(sum(SAD))
    H = 0.0
    for p in SAD:
        if p < 1.0: 
            print 'p < 1.0', p
            sys.exit()
            
        H += -(p/N)*np.log(p/N) 
    H = (np.exp(H) - 1)/(S - 1)
    return H


def simpsons_dom(SAD):
    D = 0.0
    N = sum(SAD)
    S = len(SAD)
    
    for x in SAD:
        D += x*(x-1)
    D = 1 - (D/(N*(N-1)))
    
    return D
    

def simpsons_evenness(SAD):
    D = 0.0
    N = sum(SAD)
    S = len(SAD)
    
    for x in SAD:
        D += (x*x) / (N*N)
    
    E = (1/D)/S
    
    return E


def berger_parker(SAD):
    bp = float(max(SAD))/sum(SAD)
    return bp

    
def EQ_evenness(SAD):
    
    SAD.sort()
    SAD.reverse()
    
    S = len(SAD)
    y_list = list(np.log(SAD))
    x_list = []
    for rank in range(1,S+1):
        x_list.append(rank/S)
    slope, intercept, rval, pval, std_err = stats.linregress(x_list, y_list)
    
    Eq = -2/np.pi*np.arctan(slope)
    return Eq


def e_var(SAD):
    P = np.log(SAD)
    S = len(SAD)
    X = 0
    for x in P:
        X += (x - np.mean(P))**2/S
    evar = 1 - 2/np.pi*np.arctan(X) 
    return(evar)


def get_skews(_list):
    skews = []
    for i in _list:
        skews.append(stats.skew(i))
    
    return skews


def get_modal(_list):
    
    """ Finds the mode from a kernel density function across a sample """
    exp_mode = 0.0
    density = gaussian_kde(_list)
    n = len(_list)
    xs = np.linspace(min(_list),max(_list),n)
    density.covariance_factor = lambda : .001
    density._compute_covariance()
    D = [xs,density(xs)]
    d = 0
    maxd = 0.0
    while d < len(D[1]):
        if D[1][d] > maxd:
            maxd = D[1][d]
            exp_mode = D[0][d]
        d += 1
    return exp_mode
    
def get_kdens_choose_kernel(_list,kernel):
    """ Finds the kernel density function across a sample of SADs """
    density = gaussian_kde(_list)
    n = len(_list)
    xs = np.linspace(min(_list),max(_list),n)
    #xs = np.linspace(0.0,1.0,n)
    density.covariance_factor = lambda : kernel
    density._compute_covariance()
    D = [xs,density(xs)]
    return D
    
def get_kdens(_list):
    """ Finds the kernel density function across a sample of SADs """
    density = gaussian_kde(_list)
    n = len(_list)
    #xs = np.linspace(min(_list),max(_list),n)
    xs = np.linspace(0.0,1.0,n)
    density.covariance_factor = lambda : 0.5
    density._compute_covariance()
    D = [xs,density(xs)]
    return D





def RandCompFast(Q, N, sample_size):
    
    comps = []
    while len(comps) < sample_size:

        composition = []    
        indices = []
        
        while len(indices) < N-1:
            index = random.randint(1, Q-1)
            if index in indices: continue
            else: indices.append(index)
    
        indices.sort()
        indices.append(Q)
    
        nsum = 0
        composition = [] 
        for i in indices:
            i -= sum(composition) 
            composition.append(i)
            nsum += i
        
        composition.sort()
    	composition.reverse()
        comps.append(composition)

    return comps
    
################################################################################

def get_SADs(path, dataset):

    DATA = path + '/' + dataset + '/' + dataset + '-data.txt'
    
    #print dataset
    
    minS = 10
    ct1, ct2 = 0, 0
    SADdict = {}
    
    with open(DATA) as f: 
        
        for d in f:
            if d.strip():
                sample = re.match(r'\A\S*',d).group()
                abundance = float(re.findall(r'\S*\S$',d)[0])
            
                if sample not in SADdict:
                    SADdict[sample] = [abundance]
            
                else: SADdict[sample].append(abundance)
        
    SADs = []
    SADlist = SADdict.items()
    
    for tup in SADlist:
            
        SAD = tup[1]
        if len(SAD) >= minS: 
            SAD.sort()
            SAD.reverse()
            SADs.append(SAD)
            
    return(SADs)
    
    
    
    
def GetSADsFromBiom_labeled(path, dataset):

    DATA = path + '/' + dataset + '-SSADdata.txt'
    SADdict = {}
    
    with open(DATA) as f: 
        
        for d in f:
            if d.strip():
                
                d = d.split()
                species = d[0]
                sample = d[1]
                abundance = float(d[2])
                
                if abundance > 0:
                    if sample not in SADdict:
                        SADdict[sample] = [[species, abundance]]
            
                    else: SADdict[sample].append([species, abundance])
        
            
    return SADdict
    
    
    

def getFS(Nlist, Slist, tool, zeros=False):
    
    if zeros is False: OUT = open('/Users/lisalocey/Desktop/RareBio/macrostates/'+tool+'.txt','w+')
    elif zeros is True: OUT = open('/Users/lisalocey/Desktop/RareBio/macrostates/'+tool+'weak.txt','w+')

    for i, N in enumerate(Nlist):
        S = Slist[i]
        
        print i, N, S
        N = int(N)
        S = int(S)
        
        """ use mete? """
        if tool == 'mete':
            FS_RADs = [mete.get_mete_rad(S, N)[0]]
        
        """ use compositions? """
        if tool == 'Comps':
            sample_size = 2
            FS_RADs = RandCompFast(N, S, sample_size)
        
        """ use partitions? """
        if tool == 'Parts':
            sample_size = 1
            if N > 10**4 or N/S > 50:
                FS_RADs = partitions.rand_partitions(N, S, sample_size, exact='yes', method='multiplicity')
            else: FS_RADs = partitions.rand_partitions(N, S, sample_size, exact='yes', method='divide_and_conquer')
            
        """ use random fraction? """
        #sample_size = 1
        #rel = False
        #FS_RADs = ModelsN.DomPreInt(N, sample_size, rel)                
        
        """ use sim-based log-normal? """
        #sample_size = 10
        #rel = False
        #FS_RADs = ModelsN.SimLogNormInt(N, sample_size, rel)                                       
        
        """ use sim-based Pareto? """
        #sample_size = 10
        #rel = False
        #FS_RADs = ModelsN.SimParetoInt(N, sample_size, rel)
        
        """ use sim-based Dominance Decay (int)?"""
        #sample_size = 10
        #rel = False
        #FS_RADs = ModelsN.DomDecayInt(N, sample_size, rel)
        
        for rad in FS_RADs:
            if sum(rad) != N or len(rad) != S:
                print sum(rad), N, len(rad), S
                print 'incorrect N and S'
                sys.exit()
                
            print>>OUT, rad
            
    OUT.close()
    print 'getFS(): done'
            