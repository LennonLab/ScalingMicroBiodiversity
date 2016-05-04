from __future__ import division
import  matplotlib.pyplot as plt

import pandas as pd
import linecache
import numpy as np
import os
import sys

#import statsmodels.stats.api as sms
#import statsmodels.api as sm
import statsmodels.formula.api as smf
#from statsmodels.sandbox.regression.predstd import wls_prediction_std
from statsmodels.stats.outliers_influence import summary_table

mydir = os.path.expanduser("~/GitHub/MicrobialScaling/")
mydir2 = os.path.expanduser("~/")

OUT = open(mydir + 'fig-data/Fig1_Data-Stratified.txt','w+')
print>>OUT,'DataSet, Kind, N, Nmax, S'



def Fig1_data(condition, ones):

    numMic = 0
    numMac = 0

    tail = str()
    if ones is False:
        tail = '-SADMetricData_NoMicrobe1s.txt'
    elif ones is True:
        tail = '-SADMetricData.txt'

    datasets = []
    GoodNames = []
    emp = str()

    if condition == 'open': emp = 'EMPopen'
    elif condition == 'closed': emp = 'EMPclosed'

    GoodNames = [emp, 'HMP', 'BIGN', 'TARA', 'BOVINE', 'HUMAN', 'LAUB', 'SED', 'CHU', 'CHINA', 'CATLIN', 'FUNGI', 'HYDRO', 'BBS', 'CBC', 'MCDB', 'GENTRY', 'FIA'] # all microbe data is MGRAST

    for name in os.listdir(mydir +'data/micro'):
        if name in GoodNames: pass
        else: continue

        path = mydir+'data/micro/'+name+'/'+name+tail
        num_lines = sum(1 for line in open(path))
        datasets.append([name, 'micro', num_lines])
        print name, num_lines


    for name in os.listdir(mydir +'data/macro'):
        if name in GoodNames: pass
        else: continue

        path = mydir+'data/macro/'+name+'/'+name+'-SADMetricData.txt'
        num_lines = sum(1 for line in open(path))
        datasets.append([name, 'macro', num_lines])
        print name, num_lines


    for dataset in datasets:

        name, kind, numlines = dataset
        lines = []
        small = ['BIGN', 'BOVINE', 'CHU', 'LAUB', 'SED']
        big = ['HUMAN', 'CHINA', 'CATLIN', 'FUNGI', 'HYDRO']

        if kind == 'macro':
            lines = np.random.choice(range(1, numlines+1), 100, replace=True)
        elif name in small:
            lines = np.random.choice(range(1, numlines+1), 20, replace=True)
        elif name in big:
            lines = np.random.choice(range(1, numlines+1), 50, replace=True)
        elif name == 'TARA':
            lines = np.random.choice(range(1, numlines+1), 50, replace=True)
        else:
            lines = np.random.choice(range(1, numlines+1), 50, replace=True)
        if kind == 'micro':
            path = mydir+'data/'+kind+'/'+name+'/'+name+tail
        else:
            path = mydir+'data/macro/'+name+'/'+name+'-SADMetricData.txt'

        for line in lines:
            data = linecache.getline(path, line)
            data = data.split()

            name, kind, N, S, Var, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = data

            N = np.log10(float(N))
            Nmax = np.log10(float(Nmax))
            S = np.log10(float(S))

            print>>OUT, name, kind, N, Nmax, S

            if kind == 'micro':
                numMic += 1
            elif kind == 'macro':
                numMac += 1

    print 'microbe:',numMic
    print 'macrobe:',numMac

    return


condition = 'closed'
ones = True

Fig1_data(condition, ones)
