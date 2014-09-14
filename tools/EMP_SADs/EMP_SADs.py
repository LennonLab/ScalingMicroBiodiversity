from __future__ import division
import sys
import os

import feasible_functions as ff
import random
import numpy as np


def unlabeled():
    OUT = open('/Users/lisalocey/Desktop/RareBio/data/micro/EMP/EMP-data.txt','w+')

    RADs = ff.GetSADsFromBiom('/Users/lisalocey/Desktop/RareBio/data/micro/EMP', 'EMP')
    print 'EMP', len(RADs)

    for i, RAD in enumerate(RADs):
        
        RAD = list([x for x in RAD if x != 0])
        
        for ab in RAD:       
            print>>OUT, i, ab
            
    OUT.close()
            
            
            
def labeled():
    OUT = open('/Users/lisalocey/Desktop/RareBio/data/micro/EMP/EMP-data.txt','w+')

    RADs = ff.GetSADsFromBiom_labeled('/Users/lisalocey/Desktop/RareBio/data/micro/EMP', 'EMP')
    print 'EMP', len(RADs)

    for i, RAD in enumerate(RADs):
        
        RAD = list([x for x in RAD if x != 0])
        
        for ab in RAD:       
            print>>OUT, i, ab
            
    OUT.close()
            
