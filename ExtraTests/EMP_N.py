from __future__ import division
import  matplotlib.pyplot as plt

import numpy as np
import random
import scipy as sc
from scipy import stats

import os
import sys
from scipy.stats.distributions import t

import statsmodels.stats.api as sms
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from statsmodels.stats.outliers_influence import summary_table

import itertools as it

import pandas as pd
from math import log10
import linecache

mydir = os.path.expanduser("~/GitHub/MicrobialScaling/")
mydir2 = os.path.expanduser("~/")

#sys.path.append(mydir2 + "GitHub/DiversityTools/metrics")
#import metrics as mets


def BigN():

    datasets = []
    GoodNames = ['HMP', 'EMPclosed', 'EMPopen']

    for name in os.listdir(mydir +'data/micro'):
        if name in GoodNames: pass
        else: continue

        #path = mydir2+'data/micro/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt'
        path = mydir2+'data/micro/'+name+'/'+name+'-SADMetricData.txt'

        num_lines = sum(1 for line in open(path))
        datasets.append([name, 'micro', num_lines])
        print name, num_lines


    Nlist, Slist, Evarlist, ESimplist, klist, radDATA, BPlist, NmaxList, rareSkews, KindList, StdList = [[], [], [], [], [], [], [], [], [], [], []]

    for dataset in datasets:

        Nt = 0
        name, kind, numlines = dataset
        lines = []

        if name == 'EMPclosed' or name == 'EMPopen':
            lines = np.random.choice(range(1, numlines+1), numlines, replace=False) # 166

        elif kind == 'micro':
            lines = np.random.choice(range(1, numlines+1), numlines, replace=False) #167

        #path = mydir2+'data/'+kind+'/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt'
        path = mydir2+'data/'+kind+'/'+name+'/'+name+'-SADMetricData.txt'

        for line in lines:
            data = linecache.getline(path, line)
            radDATA.append(data)

        for data in radDATA:

            data = data.split()
            name, kind, N, S, Var, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = data

            Nt += float(N)

        print name
        print '%.2e' % Nt

    return


BigN()
