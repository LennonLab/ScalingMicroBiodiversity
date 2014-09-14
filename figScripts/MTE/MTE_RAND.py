from __future__ import division
import numpy as np
from random import randint
import matplotlib.pyplot as plt
import scipy as sc
from scipy import stats

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
masses, AvgPartSizes, mass = [[], [], 4]

while mass <= 10**15:
    mass = mass*2 # increase mass by factor of 2
    masses.append(mass) # x vals

    parts = []
    for j in range(2000): parts.append(mass ** (randint(0, 9.9999*10**4)/10**5))
    AvgPartSizes.append(np.mean(parts))
    

slope, yint, r, p, s = stats.linregress(np.log(masses), np.log(AvgPartSizes))
txt = 'slope = '+str(round(slope, 2))+', r-square = '+str(round(r**2, 3))
plt.scatter(masses, AvgPartSizes, color='0.1', alpha = 0.8, label=txt)

plt.xlim(1, 10**16), plt.ylim(1, 10**16), plt.xlabel('size', fontsize=16)
plt.ylabel('f(size)', fontsize=16)

leg = plt.legend(loc=2,prop={'size':12})
leg.draw_frame(False)


plt.xscale('log'), plt.yscale('log'), ax.tick_params(axis='both', labelsize=8)    
plt.savefig('/Users/lisalocey/Desktop/MTEvsFS.png', dpi=600, pad_inches = 0.1)
plt.show()