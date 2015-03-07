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
from statsmodels.regression.quantile_regression import QuantReg

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
NofEMP = 0
SofEMP = 0

for i, line in enumerate(lines):

    data = linecache.getline(IN, line)
    SAD = eval(data)
    S = len(SAD)
    N = int(sum(SAD))

    if S < 2 or N < 10: continue

    Slist.append(float(np.log10(S)))
    Nlist.append(float(np.log10(N)))

d = pd.DataFrame({'N': list(Nlist)})
d['S'] = list(Slist)

mod = smf.quantreg('S ~ N', d) # quantile regression
res = mod.fit(q = 0.5) # get results for 0.5 quantile
print(res.summary()) # results from quantile regression


quantiles = np.arange(.05, .96, .1)
def fit_model(q):
    res = mod.fit(q=q)
    return [q, res.params['Intercept'], res.params['N']] + \
            res.conf_int().ix['N'].tolist()

models = [fit_model(x) for x in quantiles]
models = pd.DataFrame(models, columns=['q', 'a', 'b', 'lb', 'ub'])


ols = smf.ols('S ~ N', d).fit() # linear regression on log(S) and log(N)
ols_ci = ols.conf_int().ix['N'].tolist()
olsd = dict(a = ols.params['Intercept'],
           b = ols.params['N'],
           lb = ols_ci[0],
           ub = ols_ci[1])

print ols.summary()


x = np.arange(float(d.N.min()), float(d.N.max()), 0.1)
get_y = lambda a, b: a + b * x

fig, ax = plt.subplots(figsize=(8, 6))

#for i in range(models.shape[0]):
#    y = get_y(models.a[i], models.b[i])
#    ax.plot(x, y, linestyle='dotted', color='grey')

y = get_y(olsd['a'], olsd['b'])
ax.plot(x, y, color='red', label='OLS')
ax.scatter(d.N, d.S, alpha=.2)
ax.set_xlabel('log(N)', fontsize=16)
ax.set_ylabel('log(S)', fontsize=16)
plt.show()
