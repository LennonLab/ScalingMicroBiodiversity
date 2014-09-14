#!/usr/bin/env python

from __future__ import division
import sys
import numpy as np
from scipy import stats
import random, decimal
from random import choice
import re
import math
import itertools


""" Functions for generating random integer partitions of a total q having 
    exactly n parts. """


def conjugate(partition):
    """
    Find the conjugate of an integer partition. Recoded (on 24-Apr-2013) from
    the Sage source code: www.sagenb.org/src/combinat/partition.py """
    
    if partition == []:
        return []
    else:
        l = len(partition)
        conj =  [l] * partition[-1]
        for i in xrange(l - 1, 0, -1):
            conj.extend([i] * (partition[i - 1] - partition[i]))
        return conj


def NrParts(*args):
    """ Find the number of partition for a given total q and number of parts n. Recoded
        (on 24-Apr-2013) and modified from GAP source code: www.gap-system.org
        
        Modifications for speed based on the proposition that the number of partitions of
        q having n parts is equal to the number of partitions of q-n, if n > q/2 (for odd q)
        or if n >= q/2 (for even q)
        
        Using Sage:
            number_of_partitions(1, 1) = 1 = number_of_partitions(1-1), check
            number_of_partitions(2, 1) = 1 = number_of_partitions(2-1), check
            number_of_partitions(3, 1) = 1 != number_of_partitions(3-1) = 2 , check
            number_of_partitions(4, 2) = 2 = number_of_partitions(4-2), check
            number_of_partitions(5, 2) = 2 != number_of_partitions(5-2) = 3, check
            number_of_partitions(5, 3) = 2 = number_of_partitions(5-3) = 2, check
          
            number_of_partitions(10000, 4999) = 172893842978518096747262197354965370495940268765365154597697652547323374958 
            != number_of_partitions(10000 - 4999) = 172893842978518096747262197354965370495940268765365154597697652547323374960, check
                
            number_of_partitions(10000, 5000)  = 169820168825442121851975101689306431361757683049829233322203824652329144349
            = number_of_partitions(10000-5000) = 169820168825442121851975101689306431361757683049829233322203824652329144349, check
                
            number_of_partitions(10000, 9990) = 42 = number_of_partitions(10000 - 9990), check  """
    
    if len(args) == 2: # Can we use the proposition?
        
        q = args[0]
        n = args[1]
        
        if n >= q/2.0:
            args = [q-n]
    
    numparts = 0  
    if len(args) == 1: # if we're finding p(q)
        
        q = args[0]
        
        numparts = 1                             
        p = [1]*(q+1)
        for i in range(1,q+1):
            numparts = 0
            k = 1
            l = 1                         
            while 0 <= i-(l+k):
                numparts = numparts - (-1)**k * (p[i-l] + p[i-(l+k)])
                k = k + 1
                l = l + 3*k - 2
            
            if 0 <= i-l:
                numparts = numparts - (-1)**k * p[i-l]
            p[i] = numparts
    
    elif len(args) == 2: # if we're finding p(q, n)   
        q = args[0]
        n = args[1]
        
        numparts=0
        if q == n or n == 1:
            numparts = 1
        elif q < n or n == 0:
            numparts = 0
        else:
            q1 = int(q)
            k1 = int(n)
            p = [1]*q1
        
            for i in range(2,k1+1):  
                for m  in range(i+1,q1-i+1+1):
                    p[m] = p[m] + p[m-i]
            
            numparts = p[q1-k1+1]
    
    return numparts;
    

def P(D, q, k):
    """ A function to return the number of partitions of q with k or less parts.
    
    Note:
        1. Theorem: The number of partitions of q with k or less parts equals the
           number of partitions of q with k or less as the largest part (see Bona 2006).
           This is a mathematical symmetry, i.e. congruency.
           
           Bona, M. (2006). A Walk Through Combinatorics: An Introduction to Enumeration
             and Graph Theory. 2nd Ed. World Scientific Publishing Co. Singapore.
        
        2. Proposition: The number of partitions of q with k or less parts equals the
           number of partitions of q+k with k as the largest part when k>0, i.e. P(q + k, k).
           No source, but it can be shown when enumerating the entire feasible set
           or using Sage:
                      
           Checks using Sage:
               number_of_partitions(1+0, 0) = 0 = Partitions(1, max_part=0).cardinality()
               number_of_partitions(1+1, 1) = 1 = Partitions(1, max_part=1).cardinality()
               number_of_partitions(55+3, 3) = 280 = Partitions(55, max_part=3).cardinality()
               number_of_partitions(99+2, 2) = 50 = Partitions(99, max_part=2).cardinality() 
           
    Arguments:
        D : a dictionary for the number of partitions of q having k or less
            parts (or k or less as the largest part), i.e. P(q, q + k).   
        q : the total (i.e. sum across all k or n parts)
        k : the number of parts and also the size of the largest part (congruency)     
    
    Because of the congruency above (Bona 2006) k and n are sometimes used synomously.
    
    """
    
    if (q, k) not in D:
        D[(q, k)] = NrParts(q + k, k)
        # Note from the above: p(q) = p(q + q, q)
        # so NrParts(q) returns the same value as NrParts(q + q, q)
    return [D, D[(q, k)]] # return the updated dictionary and P(q + k, k).


def rand_partitions(q, n, sample_size, method='best', D={}, zeros=False):
    """
    Generate uniform random partitions of Q having N parts.
    
    Arguments:
        Q : Total sum across parts
        N : Number of parts to sum over
        sample_size : number of random partitions to generate
        method : method to use for generating the partition, options include:
            'bottom_up', 'top_down', 'divide_and_conquer', 'multiplicity', and
            'best'. Defaults to 'best'
        D : a dictionary for the number of partitions of Q having N or less
            parts (or N or less as the largest part), i.e. P(Q, Q + N). Defaults
            to a blank dictionary.
        zeros : boolean if True partitions can have zero values, if False
            partitions have only positive values, defaults to False
    
    Returns: A list of lists
    
    Notes:
        method == 'best' attempts to use the values of Q and N to infer what the 
        fastest method to compute the partition.
    
    """
    parts = []
    if zeros:
    
        """ if zeros are allowed, then we must ask whether Q >= N. if not, then
            the total Q is partitioned among a greater number of parts than there
            are, say, individuals. In which case, some parts must be zero. A random
            partition would then be any random partition of Q with zeros appended
            at the end. But, if Q >= N, then Q is partitioned among less number of
            parts than there are individuals. In which case, a random partition
            would be any random partition of Q having N or less parts. """  
        
        if q >= n:
            # if Q >= N and zero's are allowed, the first part must be equal to or less than N
            Plist = P(D, q, n) 
        elif q < n:
            # if Q < N and zero's are allowed, the first part must be equal to or less than Q
            Plist = P(D, q, q)    
    
    else:
        Plist = P(D, q - n, n) # if zero's are not allowed, the first part must be N.
        
    D = Plist[0]
    numparts = Plist[1]        
    while len(parts) < sample_size:
        rand_int = random.randrange(1, numparts + 1)
        
        if zeros:
            q1 = int(q)
            part = []
        else:
            q1 = int(q - n)
            part = [n]
        
        if method == 'bottom_up':
            part = bottom_up(part, q1, D, rand_int)
        if method == 'top_down':
            part = top_down(part, q1, D, rand_int)
        if method == 'divide_and_conquer':
            part = divide_and_conquer(part, q1, n, D, rand_int)
        if method == 'multiplicity':
            part = multiplicity(part, q1, D, rand_int)
        if method == 'best':
            if Q < 250 or N >= Q / 1.5:
                part = bottom_up(part, q1, D, rand_int)
            else:
                part = divide_and_conquer(part, q1, n, D, rand_int)
        if zeros:
            Zs = [0] * (n - len(part))
            part.extend(Zs)
        parts.append(part)
    return parts


def bottom_up(part, q, D, rand_int):
    """
    Bottom up method of generating uniform random partitions of q having n parts.
    
    Arguments:
        part : a list to hold the partition
        q : the total sum of the partition
        k : size of the largest (and also first) part
        D : a dictionary for the number of partitions of q having n or less
            parts (or n or less as the largest part), i.e. P(q + n, n).        
        rand_int : a number representing a member of the feasible set

    """    
    
    while q > 0:
        for k in range(1, q + 1): # loop through all possible values of the first/largest part
            Plist = P(D, q, k) # number of partitions of q having k or less as the largest part
            D = Plist[0]
            count = Plist[1]
            if count >= rand_int:
                Plist = P(D, q, k - 1)
                D = Plist[0]
                count = Plist[1]
                break
        part.append(k)
        q -= k
        if q == 0:
            break
        rand_int -= count
    part = conjugate(part)    
    return(part)


def top_down(part, q, D, rand_int):
    """
    Top down method of generating uniform random partitions of q having n parts.
    
    Arguments:
        part : a list to hold the partition
        q : the total sum of the partition
        k : size of the largest (and also first) part
        D : a dictionary for the number of partitions of q having n or less
            parts (or n or less as the largest part), i.e. P(q + n, n).        
        rand_int : a number representing a member of the feasible set

    """    
    
    while q > 0:
        if part != []: 
            x = min(part)
        else: 
            x = q
        for k in reversed(range(1, x + 1)): # loop through all possible values of the
        # first/largest part
            Plist = P(D, q, k) # number of partitions of q having k or less as the
            # largest part
            D = Plist[0]
            count = Plist[1]
            if count < rand_int:
                k += 1
                break
        rand_int -= count
        part.append(k)
        q -= k
    part = conjugate(part)
    return(part)


def divide_and_conquer(part, q, n, D, rand_int):
    """
    Divide and conquer method of generating uniform random partitions of q
    having n parts.
        
    Arguments:
        part : a list to hold the partition
        q : the total sum of the partition
        k : size of the largest (and also first) part
        n : number of parts to sum over
        D : a dictionary for the number of partitions of q having n or less
            parts (or n or less as the largest part), i.e. P(q + n, n).        
        rand_int : a number representing a member of the feasible set

    """
    if n >= 1 and isinstance(n, int): pass 
    else: print 'n must be a positive integer'
    
    max_int = int(n)
    min_int = int(1)
    while q > 0:
        k = random.randrange(min_int, max_int + 1) # choose a value of the 
        # largest part at random
        Plist = P(D, q, k)
        D = Plist[0]
        upper = Plist[1]
        Plist = P(D, q, k - 1)
        D = Plist[0]
        lower = Plist[1]
        if lower < rand_int and rand_int <= upper: 
            part.append(k)
            q -= k
            max_int = k
            min_int = 1
            num = int(upper - lower)
            rand_int = random.randrange(1, num + 1)
        elif rand_int > upper:
            min_int = k + 1    
        elif rand_int <= lower:
            max_int = k - 1    
    part = conjugate(part)
    return part


def get_multiplicity(q, k, D, rand_int, count): 
    """ 
    Find the number of times a value k occurs in a partition that is being
    generated at random by the multiplicity() function. The resulting
    multiplicity is then passed back to the multiplicity() function along with
    an updated value of count and an updated dictionary D
    
    Arguments:
        q : the total sum of the partition
        k : size of the largest (and also first) part 
        D : a dictionary for the number of partitions of q-k*f having k-1 or less
            parts (or k-1 or less as the largest part).                
        rand_int : a number representing a member of the feasible set
        count : number of partitions of q-k*f having k-1 or less parts 
        f : number of times k occur
        multi : list of f values of k """
        
    multi = [] # the multiplicity 
    f = 1
    while f > 0:
        Plist = P(D, (q - k * f), k - 1)
        D = Plist[0]
        count += Plist[1]
        if count >= rand_int:
            count -= Plist[1]
            multi = [k] * f
            break                
        f += 1
    return [D, count, multi]


def multiplicity(part, q, D, rand_int):
    """
    multiplicity method of generating uniform random partitions of q having n
    parts.
    
    Arguments:
        part : a list to hold the partition
        q : the total sum of the partition
        k : size of the largest (and also first) part 
        D : a dictionary for the number of partitions of q having n or less
            parts (or n or less as the largest part), i.e. P(q + n, n).        
        rand_int : a number representing a member of the feasible set

    """
    while q > 0:
        multi = []
        if part != []:
            x = min(part)
        else: 
            x = int(q)
        for k in reversed(range(1, x + 1)): # start with largest k
            Plist = P(D, q, k) # number of partitions of q having k or less as the largest part
            D = Plist[0]
            count = Plist[1]
            if count == rand_int and rand_int == 1:
                multi = [1] * q
                q = 0
                break
            if count < rand_int: # k has been found
                k += 1
                Mlist = get_multiplicity(q, k, D, rand_int, count) # now, find how many times k
                # occurs, i.e. the multiplicity of k, i.e. get_multiplicity()
                D = Mlist[0]
                count = Mlist[1]
                multi = Mlist[2]
                break
        q -= sum(multi)
        part.extend(multi)
        rand_int -= count    
    part = conjugate(part)
    return part
    
    
def test_qnk(q, n=False, k=1):
    """ A function to ensure that q, n, and k are positive integers.
    q : the total sum of the partition
    k : size of the largest (and also first) part
    n : number of parts to sum over 
    
    allows one, two, or three parameters
    by default n=q and k=1 """
    if not n: n = q
    
    if q >= 0 and isinstance(q, int): pass
    else:
        print q,'q must be a non-negative integer'
        sys.exit()
    if n >=0 and isinstance(n, int): pass
    else:
        print 'n must be a non-negative integer'
        sys.exit()
    if k >= 1 and isinstance(k, int): pass
    else:
        print k,'k must be a non-negative integer'
        sys.exit()
    
    if q < n:
        print q, n, 'q must be >= n'
        sys.exit()
    if q < k and q > 0:
        print 'q must be >= k'
        sys.exit()
    
    return
        
    
def last_lexical(q, n):
    """ Find the last lexical (i.e. most even) partition of q having n parts
    q : the total sum of the partition
    n : number of parts in the partition
    j : place holder """
    
    test_qnk(q, n)
    partition = [int(math.floor(float(q) / float(n)))] * n # a list of n integer 
    # all with the same value, but not necessarily summing to q
    _remainder = int(q % n)
    j = 0
    while _remainder > 0: # distribute the remainder evenly among parts
    # of the partition
        partition[j] += 1
        _remainder -= 1
        j += 1
    return partition


def min_max(q, n):
    """ Find the smallest possible maximum value for the first part in a partition
    of q having n parts """
    
    test_qnk(q, n)
    min_int = int(math.floor(float(q) / float(n)))
    if int(q % n) > 0:
        min_int +=1
    return min_int

    
def first_lexical(q, n, k):
    """ Find the first lexical partition of q having n parts with k as the largest part
    q : the total sum of the partition
    k : size of the largest (and also first) part 
    n : number of parts in the partition """
    
    test_qnk(q, n)
    partition = []
    if k == None:
        partition.append(q - n + 1)
        ones = [1] * (n - 1)
        partition.extend(ones)
        return partition
    
    elif k < min_max(q, n):
        return None
        
    else:
        partition.append(k)
        q -= k
        n -= 1
        while q > 0:
            k = min(k, q - n + 1)
            partition.append(k)
            q -= k
            n -= 1
        
    return partition


def next_restricted_part(partition):
    """ Find the next lexical partition of q having n parts
    q : the total sum of the partition
    n : number of parts in the partition
    halves : two halves of the partition """
    
    q = sum(partition)
    n = len(partition)
    test_qnk(q, n)
    
    if partition == last_lexical(q, n):
        return first_lexical(q, n, None)

    for i in enumerate(reversed(partition)):
        if i[1] - partition[-1] > 1:
            if i[0] == (n - 1):
                partition = first_lexical(q, n, int(i[1] - 1))
                return partition
            else:
                halves = np.split(partition, [n - i[0] - 1])
                h1 = list(halves[0])
                h2 = list(halves[1])
                next = list(first_lexical(int(sum(h2)), int(len(h2)), int(h2[0]) - 1))
                return h1 + next
