import  matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import sys

import statsmodels.stats.api as sms
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from statsmodels.stats.outliers_influence import summary_table


def func(row):
    ssad = row[0]
    ssad = ssad.split()
    ssad = map(int, ssad)

    for i, val in enumerate(ssad):
        mylist[i].append(val)

mydir = os.path.expanduser('~/GitHub/MicrobialScaling/data/micro/TARA/')
dat = pd.read_csv(mydir + '/TARA-table.txt')

names = list(dat)[0].split()
print len(names), names[0]

mylist = [list([]) for _ in xrange(len(names))]
dat.apply(func, axis=1, raw=True)

print len(mylist)
print mylist[0][0:3]

OUT = open(mydir + 'TARA-data.txt','w+')

for i, sad in enumerate(mylist):
    sad = list([int(x) for x in sad if x > 0])
    sad.sort()
    sad.reverse()
    for j, ab in enumerate(sad):
        print>>OUT, names[i], j, ab
