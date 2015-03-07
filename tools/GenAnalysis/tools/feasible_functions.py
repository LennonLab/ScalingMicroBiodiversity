from __future__ import division
import sys
import numpy as np
from scipy import stats
from scipy.stats import gaussian_kde
import re

########################################################################################################
######   A Section devoted to evenness indices and descriptive statistical functions ###################
minS = 7

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
        """Return computed Gini coefficient of inequality. This function was found at http://econpy.googlecode.com/svn/trunk/pytrix/utilities.py """

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
    if E > 1.0:
        print 'Simpsons violation',E
        print N,S, SAD
        sys.exit()
    
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

 
#######################################################################################################

def get_SADs(path, dataset):
    
    minS = 10
    DATA = open(path + '/' + dataset + '.txt','r')
    ct1 = 0
    ct2 = 0
    d = DATA.readline()
    m0 = re.match(r'\A\S*',d).group()
    site_name = str(m0)
    m2 = float(re.findall(r'\S*\S$',d)[0])
    
    SAD = [float(m2)]
    SADs = []
        
    for d in DATA:
        ct1+=1
        m1 = re.match(r'\A\S*',d).group()
        if m1 == m0:
            m2 = float(re.findall(r'\S*\S$',d)[0])
            if m2 > 0:
                SAD.append(m2)
                
        else:
            site_name = m0
            m0 = m1
            if len(SAD) >= minS:
                SAD.sort()
                SAD.reverse()
                SADs.append([site_name, SAD]) # can also append, site_name, len(SAD), and sum(SAD)
                
                ct2+=1
            SAD = []
            abundance = re.findall(r'\S*\S$',d)[0]
            
            if abundance > 0:SAD.append(float(abundance))
    
    #SAD.sort()
    #SAD.reverse()
    #SADs.append(SAD)
    
    DATA.close()
    SADs.append([site_name, SAD])
    return(SADs)
    
    

def get_hottest_SAD(unique_SADs):
    """ Find the SAD in a random sample with the greatest average commonness 
        among its ranked abundance states. This SAD is taken to represent the 
        central tendency of the set, based on the SAD shape. """
    
    #if len(unique_SADs) > 500:
        #unique_SADs = random.sample(unique_SADs,500)
    
    N = sum(unique_SADs[0])
    S = len(unique_SADs[0])
    a1 = 0 # SAD mean
    v1 = 0 # SAD variance
    for rad in unique_SADs:
        in_common = []
        ct1 = 0
        for a in rad: # for each rank
            c = 0
            for sad in unique_SADs: 
                if a == sad[ct1]:
                    c += 1
            in_common.append(np.log(c))
            ct1 += 1
        a2 = np.mean(in_common)
        v2 = np.var(in_common, ddof=1)  
        if a2 > a1:
            a1 = a2
            v1 = v2
            xRAD = rad
        elif a2 == a1:
            if v2 < v1:
                a1 = a2
                v1 = v2
                xRAD = rad
    #percentile_evar = stats.percentileofscore(sample_evar,obs_evar)
    return xRAD
    
    