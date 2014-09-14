#!/usr/bin/env python

from __future__ import division
import sys
import os
import partitions as parts
import time

""" This file generates .txt files holding times taken for partitioning algorithms to generate random partitions for different combinations of Q and N.
    The time results held in these files were used to generate figure ? of Locey and McGlinn 2013. """

algorithms = ['multiplicity','top_down','divide_and_conquer','bottom_up']

sample_size = 300 # takes approx 300 random macrostates to capture shape of the feasible set

qs = range(20,100,5)
qs2 = range(100,1020,20)
qs.extend(qs2)
zeros = [True, False]

for zero in zeros:
    for alg in algorithms:
    
        for q in qs:
            
            if q <= 100: step = 1
            elif q <= 300: step = 5
            elif q <= 600: step = 10
            else: step = 20
            
            n = int(step)
                
            while n <= q:  
                print alg, zero, alg, q
                
                times = []
                D = {}    
                t0 = time.time()
                x = parts.rand_partitions(q, n, sample_size, alg, D, zeros)
                t = time.time() - t0
                times.append([n,t])
                
                if zero == True:
                    OUT = open('timeFiles/Python_' + alg + '_zeros_q=' + str(q) + '.txt','a+')
                elif zero == False:
                    OUT = open('timeFiles/Python_' + alg + '_q=' + str(q) + '.txt','a+')
                
                for i in times:
                    print>>OUT, i[0], i[1]
                
                OUT.close()
                n+=step
                
