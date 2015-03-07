from __future__ import division

import sys
import macroecotools
import macroeco_distributions
import mete
import cloud

import  matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import random



def RandomComposition(Q,N):
    composition = []
    
    indices = random.sample(range(1,Q+1),N-1) # get a random sample of k-1 objects from n-1 objects
    indices.sort()
    Qlist = [1]*Q
    
    partsOF = [Qlist[i:j] for i, j in zip([0]+indices, indices+[None])]
    for p in partsOF:
        composition.append(sum(p))
    
    return composition
    
    
def get_LowerTrunc_GeomSeries(N,S):
    
    rank = range(1,S+1)
    cdf = [(S-i+0.5)/S for i in rank]
    SNratio = S/N
    abd = macroeco_distributions.trunc_geom.ppf(np.array(cdf), SNratio, N)
    
    return abd


def geoS(N, S):
    job_id = cloud.call(get_LowerTrunc_GeomSeries, N, S, _type='m1')
    geomSeries = cloud.result(job_id)   
    geomSeries = get_LowerTrunc_GeomSeries(N, S)
    geomSeries = np.log(geomSeries)
    plt.plot(geomSeries, color='m',lw=3,label='Geometric series\nN='+str(N)+' S='+str(S))
    print 'geometric series: done'
    return

def logS(N, S):
    job_id = cloud.call(mete.get_mete_rad, S, N, _type='m1')
    logSeries = cloud.result(job_id)   
    #logSeries = mete.get_mete_rad(S, N) # The expected SAD from the random sample
    logSeries = np.log(logSeries[0])
    plt.plot(logSeries, color='gray', lw=3, label='Log-series\nN='+str(N)+' S='+str(S))
    print 'log-series: done'
    return

def randGeoS(N, S, n):
    
    for n in range(n):
        #job_id = cloud.call(RandomComposition, N, S, _type='m1')
        #RAD = cloud.result(job_id)
        RAD = RandomComposition(N, S)
        RAD.sort()
        RAD.reverse()
        RAD = np.log(RAD)
        plt.plot(RAD, color='gray', lw=2, label='Random Geometric\nN='+str(N)+' S='+str(S))
    print 'Compositions: done'
    return


fig = plt.figure()
ax = fig.add_subplot(1,1,1)

N = 10**20
S = 10**3
n = 10

randGeoS(N, S, n)        


plt.show()