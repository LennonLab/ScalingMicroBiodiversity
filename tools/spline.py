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


OrC = 'closed'
path = mydir2 +'data/micro/EMP' + OrC
IN = path + '/EMP' + OrC + '-SADs.txt'

num_lines = sum(1 for line in open(IN))
lines = sample(range(1, num_lines+1), 1000)

S = []
N = []

for i, line in enumerate(lines):

    d = linecache.getline(IN, line)
    SAD = eval(d)
    s = len(SAD)
    if s < 1: continue
    n = int(sum(SAD))

    S.append(float(np.log(s)))
    N.append(float(np.log(n)))


data = pd.DataFrame({"N": N})
data["S"] = list(S)

plt.title("B-spline basis example (degree=3)");

#N = np.linspace(1., 10000000., 100)
S = patsy.dmatrix("S ~ bs(N, df=3)", {"N": N, "S": S})

# Define some coefficients
b = np.array([1.3, 0.6, 0.9])

# Plot B-spline basis functions (colored curves) each multiplied by its coeff
plt.plot(N, y)

# Plot the spline itself (sum of the basis functions, thick black curve)
plt.plot(N, np.dot(y, b), color='k', linewidth=3)
plt.show()
