from __future__ import division
#from __future__ import print_function

import sys
import os
import linecache

from random import sample

import patsy
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
#from statsmodels.regression.quantile_regression import QuantReg

mydir = os.path.expanduser("~/Desktop/Repos/rare-bio/")
sys.path.append(mydir)
mydir2 = os.path.expanduser("~/Desktop/")
sys.path.append(mydir2)


OrC = 'open'
path = mydir2 +'data/micro/EMP' + OrC
IN = path + '/EMP' + OrC + '-SADs.txt'

num_lines = sum(1 for line in open(IN))
lines = sample(range(1, num_lines+1), 1000)

Slist = []
Nlist = []

for i, line in enumerate(lines):

    data = linecache.getline(IN, line)
    SAD = eval(data)
    S = len(SAD)
    if S < 1: continue
    N = int(sum(SAD))

    Slist.append(float(np.log(S)))
    Nlist.append(float(np.log(N)))

d = pd.DataFrame({'N': list(Nlist)})
d['S'] = list(Slist)


mod1 = smf.quantreg('S ~ N', d)
mod2 = smf.quantreg('S ~ N + I(N ** 2.0)', d)

quantiles = np.array([.99])

def fit_model(q, degree):

    if degree == 1: res = mod1.fit(q=q)
    elif degree == 2: res = mod2.fit(q=q)
    return res



fig = plt.figure(figsize=(8, 6))

for i in range(2, 3):

    fig.add_subplot(1, 1, 1)

    quant = fit_model(0.95, i)
    print quant.summary()
    qy = quant.fittedvalues.tolist()

    if i == 1: ols1 = smf.ols('S ~ N', d).fit()
    elif i == 2: ols1 = smf.ols('S ~ N + I(N ** 2.0)', d).fit()

    y = ols1.fittedvalues.tolist()

    plt.scatter(Nlist, y, color='red', label='OLS')
    plt.scatter(Nlist, qy, color='grey', label='99th')
    plt.scatter(Nlist, Slist, alpha=.2)

    plt.xlabel('log(N)', fontsize=16)
    plt.ylabel('log(S)', fontsize=16)
    plt.legend(loc=2, fontsize=10)


plt.show()
