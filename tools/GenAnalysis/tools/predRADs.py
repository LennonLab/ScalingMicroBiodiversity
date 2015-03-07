from __future__ import division
import sys

sys.path.append("/Users/lisalocey/Desktop/global/GenAnalysis/tools/")
import macroeco_distributions as md
import macroecotools as mt
import feasible_functions as ff
import mete
import pln

import cloud
import numpy as np
import random
import math


def truncate(RAD, maxn):

    while max(RAD) > maxn:
        m = max(RAD)
        i = RAD.index(m)
        n1 = RAD[i]
        n2 = random.randint(1,n1-1)
        RAD[i] -= n2

        i = random.randint(0, len(RAD)-1)
        RAD[i] += n2

    return RAD


def get_GeomSeries(N,S,zeros=False):

    rank = range(1,S+1)
    cdf = [(S-i+0.5)/S for i in rank]
    SNratio = S/N
    if zeros == False:
        abd = md.trunc_geom.ppf(np.array(cdf), SNratio, N)

    return abd


def RandComp(Q, N, maxn=False):

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

    RAD = truncate(composition, maxn)
    RAD.sort()
    RAD.reverse()

    return composition



#def PowerFrac(N, S, maxn):



def RandFracPre(N, S, maxn):

    RAD = [N]
    while len(RAD) < S:

        i = random.randint(0, len(RAD)-1)
        n1 = RAD[i]
        if n1 < 2: continue
        n2 = random.randint(1,n1-1)
        RAD[i] -= n2
        RAD.append(n2)

    RAD = truncate(RAD, maxn)
    RAD.sort()
    RAD.reverse()

    return RAD

def RandFrac(N, S, maxn):

    RAD = [N]
    while len(RAD) < S:

        i = random.randint(0, len(RAD)-1)
        n1 = RAD[i]
        if n1 < 2: continue
        n2 = random.randint(1,n1-1)
        RAD[i] -= n2
        RAD.append(n2)

    RAD = truncate(RAD, maxn)
    RAD.sort()
    RAD.reverse()

    return RAD


def getPred(N, S, maxn, model, sample_size=1):

    RADs = []
    for i in range(sample_size):

        if model == 'compositions':
            # Predicted from compositions (basically geometric series)
            predRAD = predRADs.RandComp(RAD)

        elif model == 'power fraction':
            # Predicted from Fraction 1: Power Fraction
            predRAD = predRADs.PowerFrac(RAD)

        elif model == 'random fraction non-preemption':
            # Predicted from Fraction 2: Random non-preemption
            predRAD = predRADs.RandFracPre(N, S, 'pln')

        elif model == 'random fraction':
            # Predicted from Fraction 3: Random preemption
            predRAD = predRADs.RandFrac(N, S)

        RADs.append(predRAD)

    predRAD = ff.get_hottest_SAD(RADs)

    return predRAD
