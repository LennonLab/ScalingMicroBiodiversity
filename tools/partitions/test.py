#! /usr/bin/env python

from __future__ import division
import sys
import os
import partitions as parts
sys.path.append("/home/kenlocey/partitions/new/metrics")
import metrics as mt
from os import path, access, R_OK  # W_OK for write permission
from scipy.stats import gaussian_kde
from pylab import *
import numpy as np
from scipy import stats
from scipy.stats import kstest
import random
from random import choice
import re
import math
import random, decimal


def test_first_lexical():
    
    print '\nTesting to ensure the first lexical partition is correctly generated.'
    partition = [123]
    q = 123
    n = 1
    k = 123
    first_partition = parts.first_lexical(q, n, k)
    if partition != first_partition:
        print 'first_lexical() is broken. Test 1 FAIL'
    else:
        print 'first_lexical() works. Test 1 PASS'
    
    partition = [1,1,1,1,1,1,1,1,1,1]
    q = sum(partition)
    n = len(partition)
    k = max(partition)
    first_partition = parts.first_lexical(q, n, k)
    if partition != first_partition:
        print 'first_lexical() is broken. Test 2 FAIL'
    else:
        print 'first_lexical() works. Test 2 PASS'
    
    partition = [10, 10, 1, 1, 1, 1]
    q = sum(partition)
    n = len(partition)
    k = max(partition)
    first_partition = parts.first_lexical(q, n, k)
    if partition != first_partition:
        print 'first_lexical() is broken. Test 3 FAIL'
    else:
        print 'first_lexical() works. Test 3 PASS'
    
    return

test_first_lexical()
    

def test_last_lexical():
    
    print '\nTesting to ensure the last lexical partition is correctly generated.'
    partition = [1]
    q = sum(partition)
    n = len(partition)
    last = parts.last_lexical(q, n)
    if partition != last:
        print 'last_lexical is broken. Test 1 FAIL'
    else:
        print 'last_lexical works. Test 1 PASS' 
    
    partition = [1,1,1,1,1,1,1]
    q = sum(partition)
    n = len(partition)
    last = parts.last_lexical(q, n)
    if partition != last:
        print 'last_lexical is broken. Test 2 FAIL'
    else:
        print 'last_lexical works. Test 2 PASS' 
    
    partition = [11,11,11,10,10,10,10]
    q = sum(partition)
    n = len(partition)
    last = parts.last_lexical(q, n)
    if partition != last:
        print 'last_lexical is broken. Test 3 FAIL'
    else:
        print 'last_lexical works. Test 3 PASS' 
        
    return

test_last_lexical()


def test_next_restricted_part():
    
    print '\nTesting to ensure that the next lexical partition is correctly generated.'
    partition = [10,9,8,1]
    answer = [10,9,7,2]
    q = sum(partition)
    n = len(partition)
    next = parts.next_restricted_part(partition)
    if next != answer:
        print 'next_restricted_part is broken. Test 1 FAIL'
    else:
        print 'next_restricted_part works. Test 1 PASS'
        
    return

test_next_restricted_part()


def test_min_max():
    
    print '\nTesting to ensure the smallest maximum part is correctly calculated.'
    q = 100
    n = 5
    minmax1 = 20
    minmax2 = parts.min_max(q, n)
    if minmax1 != minmax2:
        print 'min_max() is broken. Test 1 FAIL'
    else:
        print 'min_max() works. Test 1 PASS'
        
    q = 100
    n = 3
    minmax1 = 34
    minmax2 = parts.min_max(q, n)
    if minmax1 != minmax2:
        print 'min_max() is broken. Test 2 FAIL'
    else:
        print 'min_max() works. Test 2 PASS'
   
    return

test_min_max()


def test_P():
    
    print '\nTesting to ensure that the number of partitons of a total q having k or less as the largest part is correctly calculated.'
    D= {}
    q = 0
    k = 0
    answer = 1 # by convention
    numparts = parts.P(D, q, k)
    if numparts[1] != answer:
        print 'P() is broken. Test 1 FAIL'
    else:
        print 'P() works. Test 1 PASS'
    
    q = 0
    k = 123
    answer = 1 # by convention
    numparts = parts.P(D, q, k)
    if numparts[1] != answer:
        print 'P() is broken. Test 2 FAIL'
    else:
        print 'P() works. Test 2 PASS'
        
    q = 123
    k = 123
    answer = 2552338241
    numparts = parts.P(D, q, k)
    if numparts[1] != answer:
        print 'P() is broken. Test 3 FAIL'
    else:
        print 'P() works. Test 3 PASS'
        
    q = 123
    k = 0
    answer = 0
    numparts = parts.P(D, q, k)
    if numparts[1] != answer:
        print 'P() is broken. Test 4 FAIL'
    else:
        print 'P() works. Test 4 PASS'
        
    q = 12345
    k = 123
    answer = 488259162924433580696194373878466788895554319556195978121822180221785381227453675217103501271281020550492
    numparts = parts.P(D, q, k)
    if numparts[1] != answer:
        print 'P() is broken. Test 5 FAIL'
    else:
        print 'P() works. Test 5 PASS'
    
    return

test_P()

def test_NrParts():
    
    print '\nTesting to ensure the number of partitions for a given total (q) and number of parts (n) is correctly calculated.'
    q = 0
    numparts = parts.NrParts(q)
    if numparts != 1:
        print 'NrParts is broken. p(0) = 1 !=',numparts,'Test 1 FAIL'
    else:
        print 'NrParts works. Test 1 PASS'
    
    n = 1
    numparts = parts.NrParts(q, n)
    if numparts != 1:
        print 'NrParts is broken. p(0, 1) = 1 !=',numparts,'Test 2 FAIL'
    else:
        print 'NrParts works. Test 2 PASS'
        
    q = 12345
    answer = 69420357953926116819562977205209384460667673094671463620270321700806074195845953959951425306140971942519870679768681736
    numparts = parts.NrParts(q)
    if numparts != answer:
        print 'NrParts is broken. Test 3 FAIL'
    else:
        print 'NrParts works. Test 3 PASS'
        
    n = 876
    answer = 2513021291084594958506275164399019380305807385021575642918576266152785748154845293483022243930843412000601280303427
    numparts = parts.NrParts(q, n)
    if numparts != answer:
        print 'NrParts is broken. Test 4 FAIL'
    else:
        print 'NrParts works. Test 4 PASS'
        
    return

test_NrParts()


def test_conjugate():
    
    print '\nTesting to ensure partitions are correctly conjugated.'
    partition = [1]
    conjugate = [1]
    conj = parts.conjugate(partition)
    if conj != conjugate:
        print 'parts.conjugate() is broken. Test 1 FAIL '
    else:
        print 'parts.conjugate() works.  Test 1 PASS'
        
    partition = [10]
    conjugate = [1,1,1,1,1,1,1,1,1,1]
    conj = parts.conjugate(partition)
    if conj != conjugate:
        print 'parts.conjugate() is broken. Test 2 FAIL '
    else:
        print 'parts.conjugate() works. Test 2 PASS '
        
    partition = [3,2,1] # symmetrical partition and conjugate
    conjugate = [3,2,1]
    conj = parts.conjugate(partition)
    if conj != conjugate:
        print 'parts.conjugate() is broken. Test 3 FAIL'
    else:
        print 'parts.conjugate() works. Test 3 PASS'
    
    partition = [12,9,8,8,7,5,4,2,2,1]
    conjugate = [10, 9, 7, 7, 6, 5, 5, 4, 2, 1, 1, 1]
    conj = parts.conjugate(partition)
    if conj != conjugate:
        print 'parts.conjugate() is broken. Test 4 FAIL'
    else:
        print 'parts.conjugate() works. Test 4 PASS'
    
    return

test_conjugate()


def get_kdens(summands):
    """ Finds the kernel density function across a sample of parts
    of partitions for a given total (N) and number of parts (S) """
    
    density = gaussian_kde(summands)
    n = len(summands)
    xs = np.linspace(float(min(summands)),float(max(summands)),n)
    density.covariance_factor = lambda : .5
    density._compute_covariance()
    D = [xs,density(xs)]
    return D


def bias_check():

    print '\nTesting algorithms for bias across combinations of q and n, allowing or excluding zero-valued parts.'
    qn_combos = [[50,10],[100,20],[200,40]]
    
    for combo in qn_combos:
        q = combo[0]
        n = combo[1]
        
        sagepartitions = []     
        DATA = open('/home/kenlocey/partitions/new/testfiles/sage_zeros_q=' + str(q) + '_n='+str(n)+'.txt','r')
        for line in DATA:
            partition = eval(line)
            sagepartitions.append(partition)
        sagevars = []
        for partition in sagepartitions:
            partition = np.log(partition)
            var = np.var(partition, ddof=1)
            sagevars.append(var)
              
        names = ['divide_and_conquer','multiplicity','top_down','bottom_up']
        for name in names:
            D = {}
            sample_size = len(sagepartitions)
            
            passes = 0
            tries = 0
            while tries < 10:
                partitions = parts.rand_partitions(q, n, sample_size, name, D, zeros=True)
                #partitions = [list(x) for x in set(tuple(x) for x in partitions)]
            
                myvars = []
                for partition in partitions:
                    partition = np.log(partition)
                    var = np.var(partition, ddof=1)
                    myvars.append(var)
                
                tries+=1
                ks, p = stats.ks_2samp(sagevars, myvars)
                #t, p = stats.ttest_ind(sagevars, myvars)
                if p > 0.05:
                    passes+=1
            print name,'(q =', q,' n =',n,') with zeros: PASSED', passes,'out of',tries 
            
        DATA.close()  
        
        sagepartitions = []     
        DATA = open('/home/kenlocey/partitions/new/testfiles/sage_q=' + str(q) + '_n='+str(n)+'.txt','r')
        for line in DATA:
            partition = eval(line)
            sagepartitions.append(partition)
        sagevars = []
        for partition in sagepartitions:
            partition = np.log(partition)
            var = np.var(partition, ddof=1)
            sagevars.append(var)
              
        names = ['divide_and_conquer','multiplicity','top_down','bottom_up']
        for name in names:
            D = {}
            sample_size = len(sagepartitions)
            
            passes = 0
            tries = 0
            while tries < 10:
                partitions = parts.rand_partitions(q, n, sample_size, name, D, zeros=False)
                #partitions = [list(x) for x in set(tuple(x) for x in partitions)]
            
                myvars = []
                for partition in partitions:
                    partition = np.log(partition)
                    var = np.var(partition, ddof=1)
                    myvars.append(var)
                
                tries+=1
                ks, p = stats.ks_2samp(sagevars, myvars)
                #t, p = stats.ttest_ind(sagevars, myvars)
                if p > 0.05:
                    passes+=1
            print name,'(q =', q,' n =',n,') no zeros: PASSED', passes,'out of',tries
            
        DATA.close()           
    return
    
bias_check()


def find_all():
    
    print '\nTesting that each random partitioning algorithm can discover the entire feasible set.'
    q = 20
    n = 5
    answer = 84 # there 84 partitions of 20 having 5 parts
    sample_size = 1000
    names = ['divide_and_conquer','multiplicity','top_down','bottom_up']
    for name in names:
        passtest = 0
        feasibleset = []
        D = {}
        while passtest < 1:
            partitions = parts.rand_partitions(q, n, sample_size, name, D, zeros=False)
            partitions = [list(x) for x in set(tuple(x) for x in partitions)]
            feasibleset.extend(partitions)
            feasibleset = [list(x) for x in set(tuple(x) for x in feasibleset)]
            if len(feasibleset) == answer:
                print 'entire feasible set found using',name,'PASS'
                passtest +=1
            
            elif len(feasibleset) > answer:
                print 'feasible set is larger than it should be. FAIL.'
                break
            else:
                print 'working to randomly generate entire feasible set for q='+str(q)+' n='+str(n)
            
    return
               
find_all()  

    
print 'TESTS FINISHED'
