from __future__ import division
import numpy as np
from random import randint, randrange
import matplotlib.pyplot as plt
import scipy as sc
from scipy import stats

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
masses, mass = [[], 3]

vals = [640, 700, 760, 820, 880, 940, 980]
colors = ['m', 'k', 'b', 'r', 'g', 'c', 'orange']

while mass <= 10**15:
    mass = mass*2 # increase mass by factor of 2
    masses.append(mass) # x vals

for i, val in enumerate(vals):
    AvgPartSizes = []
    for mass in masses:
        parts = []
        for j in range(1000): parts.append(mass ** (randint(0, val)/10**3))
        AvgPartSizes.append(np.mean(parts))
    
    slope, yint, r, p, s = stats.linregress(np.log(masses), np.log(AvgPartSizes))
    txt = 'upper lim='+str(val/1000)+' slope='+str(round(slope, 2))+', r^2='+str(round(r**2, 2))

    plt.scatter(masses, AvgPartSizes, color=colors[i], alpha = 0.4, label=txt)

leg = plt.legend(loc=2,prop={'size':12})
leg.draw_frame(False)

plt.xlim(1, 10**16), plt.ylim(1, 10**16), plt.xlabel('size', fontsize=16)
plt.ylabel('f(size)', fontsize=16)
plt.xscale('log'), plt.yscale('log'), ax.tick_params(axis='both', labelsize=8)    

plt.savefig('/Users/lisalocey/Desktop/MTEcolors.png', dpi=600, pad_inches = 0.1)
plt.show()