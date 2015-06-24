from __future__ import division
import numpy as np

import os
import linecache

mydir = os.path.expanduser("~/GitHub/rare-bio/")
mydir2 = os.path.expanduser("~/")



def combos():

    NDcombos = []
    datasets = []
    BadNames = ['.DS_Store', 'BCI', 'AGSOIL', 'SLUDGE', 'NABC']
    radDATA = []

    for name in os.listdir(mydir2 +'data/micro'):
        if name in BadNames: continue

        #path = mydir2+'data/micro/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt'
        path = mydir2+'data/micro/'+name+'/'+name+'-SADMetricData.txt'

        num_lines = sum(1 for line in open(path))
        datasets.append([name, 'micro', num_lines])
        print name, num_lines

    for name in os.listdir(mydir2 +'data/macro'):
        if name in BadNames: continue

        #path = mydir2+'data/macro/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt'
        path = mydir2+'data/macro/'+name+'/'+name+'-SADMetricData.txt'

        num_lines = sum(1 for line in open(path))
        datasets.append([name, 'macro', num_lines])
        print name, num_lines

    for dataset in datasets:

        name, kind, numlines = dataset
        lines = []

        if name == 'EMPclosed' or name == 'EMPopen':
            lines = np.random.choice(range(1, numlines+1), 25, replace=True)

        if kind == 'micro': lines = np.random.choice(range(1, numlines+1), 50, replace=True)
        else: lines = np.random.choice(range(1, numlines+1), 80, replace=True)

        #path = mydir2+'data/'+kind+'/'+name+'/'+name+'-SADMetricData_NoMicrobe1s.txt'
        path = mydir2+'data/'+kind+'/'+name+'/'+name+'-SADMetricData.txt'

        for line in lines:
            data = linecache.getline(path, line)
            radDATA.append(data)

    for data in radDATA:

        data = data.split()
        name, kind, N, S, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = data

        N = int(float(N))
        S = int(float(S))

        if S < 2 or N < 10: continue

        Nmax = int(float(Nmax))

        NDcombos.append([N, Nmax])

    return NDcombos


NDcombos = combos()
print NDcombos
