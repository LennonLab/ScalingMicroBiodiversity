from __future__ import division
import  matplotlib.pyplot as plt

import sys
import os
#sys.path.append("/path")
import feasible_functions as ff
#from os import path, access, R_OK  # W_OK for write permission
import random
import numpy as np


def AgAndBCI():

    datasets = ['AGSOIL', 'BCI']
    
    colors = ['Yellow', 'Cyan']
    
    labels = ['Ag. Soils', 'Trees on BCI']

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    ax.spines['top'].set_color('w')
    ax.spines['bottom'].set_color('w')
    ax.spines['left'].set_color('w')
    ax.spines['right'].set_color('w')        
    ax.xaxis.label.set_color('w')
    ax.tick_params(axis='both', colors='w', length=8, width=4)
        
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(4)
    
    for j, name in enumerate(datasets):
    
        if j == 0:
            path = '/Users/lisalocey/Desktop/RareBio/data/micro'
            RADs = ff.get_SADs(path, name)
        
        elif j == 1:
            path = '/Users/lisalocey/Desktop/RareBio/data/macro'
            RADs = ff.get_SADs(path, name)    
        
        color = colors[j]
        
        N = 0
        S = 0
        
        print name, len(RADs)
    
        for RAD in RADs:
        
            N = sum(RAD)
            S = len(RAD)
            
            RAD.sort()
            RAD.reverse()
        
            RAD = list([x for x in RAD if x != 0])
            #RAD = list(np.log(RAD))
        
            S = len(RAD)
        
            if S >= 50:
                ranks = range(1,len(RAD)+1)
                for i, rank in enumerate(ranks):
                    ranks[i] = rank # /S
            
                if j == 1: plt.plot(ranks, RAD, color = color, lw = 3)
                else: plt.plot(ranks, RAD, color = color, lw = 2)
    
        plt.plot([0],[-1], c = color, lw = 2, label=labels[j])#+'\nN='+str(int(N))+', S='+str(int(S)))
        #plt.xlim(1, 400)
        plt.yscale('log')
        #plt.xscale('log')
        #plt.ylim(0,12.0)
        leg = plt.legend(loc=1, prop={'size':22})
        leg.draw_frame(False)
        
        leg.get_frame().set_alpha(0.0)
        for text in leg.get_texts():
            text.set_color('w')
        
        plt.ylabel('Abundance', color='w', fontsize=40)
        plt.xlabel('Rank', color='w', fontsize=40)
        plt.tick_params(axis='both', which='major', labelsize=22) 
        
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/AgAndBCI.png',
             dpi=600, transparent = True, bbox_inches = 'tight', pad_inches=0.03)
    plt.close()
    return
    









def microRADs():

    datasets = []
    for name in os.listdir('/Users/lisalocey/Desktop/RareBio/data/micro'): 
        datasets.append(name)                                 
                              
    #datasets.remove('AGSOIL')
    labels = ['Agricultural Soils','Arctic Ocean', 'Arctic Soil', 'Human Feces',
        'Indoor Fungi','Hydrothermal vents', 'Temp & Trop Soil',
        'Activated Sludge']
    
    colors = ['Yellow', 'Yellow', 'Gold', 'LemonChiffon',
                'PapayaWhip', 'Moccasin', 'LightYellow', 'LightGoldenrodYellow']
    
    #labels = datasets
    #datasets = ['AGSOIL']
    #labels = ['Agricultural soil']
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    ax.spines['top'].set_color('w')
    ax.spines['bottom'].set_color('w')
    ax.spines['left'].set_color('w')
    ax.spines['right'].set_color('w')        
    ax.xaxis.label.set_color('w')
    ax.tick_params(axis='both', colors='w', length=8, width=4)
        
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(4)
    
    ct = 0
    for j, name in enumerate(datasets):
    
        path = '/Users/lisalocey/Desktop/RareBio/data/micro'
        RADs = ff.get_SADs(path, name)
    
        #r = lambda: random.randint(0,255)
        color = colors[j] #'#%02X%02X%02X' % (r(),r(),r())
        
        N = 0
        S = 0
        
        print name, len(RADs)
    
        for RAD in RADs:
        
            N = sum(RAD)
            S = len(RAD)
            
            RAD.sort()
            RAD.reverse()
        
            RAD = list([x for x in RAD if x != 0])
            #RAD = list(np.log(RAD))
        
            S = len(RAD)
        
            if S >= 50:
                ranks = range(1,len(RAD)+1)
                for i, rank in enumerate(ranks):
                    ranks[i] = rank/S
            
                plt.plot(ranks, RAD, color = color, lw = 0.75, alpha=0.8)
                
    
        plt.plot([0],[-1], c = color, lw = 2, label=labels[j])#+'\nN='+str(int(N))+', S='+str(int(S)))
        #plt.xlim(1, 400)
        plt.yscale('log')
        #plt.xscale('log')
        #plt.ylim(0,12.0)
        leg = plt.legend(loc=1, prop={'size':18})
        leg.draw_frame(False)
        
        leg.get_frame().set_alpha(0.0)
        for text in leg.get_texts():
            text.set_color('w')
        
        #plt.ylabel('Abundance', color='w', fontsize=25)
        #plt.xlabel('Rank', color='w', fontsize=25)
        plt.tick_params(axis='both', which='major', labelsize=22) 
        
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/microbeRADs.png',
             dpi=600, transparent = True, bbox_inches = 'tight', pad_inches=0.03)
    plt.close()
    return
    
    
    
def NvS():

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    ax.spines['top'].set_color('w')
    ax.spines['bottom'].set_color('w')
    ax.spines['left'].set_color('w')
    ax.spines['right'].set_color('w')        
    ax.xaxis.label.set_color('w')
    ax.tick_params(axis='both', colors='w') 
    
    microdatasets = []
    for name in os.listdir('/Users/lisalocey/Desktop/XRADs/micro'): 
        microdatasets.append(name)                                 
        
    ct = 0
    Ns = []
    Ss = []
        
    for j, name in enumerate(microdatasets):
    
        path = '/Users/lisalocey/Desktop/XRADs/micro'
        RADs = ff.get_SADs(path, name)
    
        r = lambda: random.randint(0,255)
        color = '#%02X%02X%02X' % (r(),r(),r())
        
        for RAD in RADs:
        
            Ns.append(sum(RAD))
            Ss.append(len(RAD))

    plt.scatter(Ns, Ss, color = 'm', linewidth = 0.0, s=10, alpha=0.8)
    #plt.scatter([-1],[-1], c = 'm', s=10, label='microbes')
        
    
    macrodatasets = []
    for name in os.listdir('/Users/lisalocey/Desktop/XRADs/macro'): 
        macrodatasets.append(name)                                 
        
    ct = 0
    Ns = []
    Ss = []
    for j, name in enumerate(macrodatasets):
    
        path = '/Users/lisalocey/Desktop/XRADs/macro'
        RADs = ff.get_SADs(path, name)
    
        r = lambda: random.randint(0,255)
        color = '#%02X%02X%02X' % (r(),r(),r())
        
        for RAD in RADs:
        
            Ns.append(sum(RAD))
            Ss.append(len(RAD))

        plt.scatter(Ns, Ss, color = 'gray', linewidth = 0.0, s=10, alpha=0.8)
    #plt.scatter([-1],[-1], c = 'gray', s=10, label='macrobes')
            
    plt.yscale('log')
    plt.xscale('log')
    #leg = plt.legend(loc=2, prop={'size':16})
    #leg.draw_frame(False)
    #plt.ylabel('Richness',fontsize=25)
    #plt.xlabel('Total abundance', fontsize=25)
    plt.tick_params(axis='both', which='major', labelsize=22) 
        
    plt.savefig('/Users/lisalocey/Desktop/XRADs/eebFigs/NvS.png',
             dpi=600, transparent = False, bbox_inches = 'tight', pad_inches=0.03)
    plt.close()
    return        

        
                        
             
def macroRADs():

    datasets = []
    for name in os.listdir('/Users/lisalocey/Desktop/RareBio/data/macro'): 
        datasets.append(name)                                 
    
    labels = ['Breeding Bird', 'Tropical Trees', 'Land Birds', 
             'Forest transects', 'Forest plots', 'Mammals', 'Butterflies']
    #labels = ['Trees on Barro Colorado Island']
    #labels = datasets
    
    colors = ['Cyan', 'Aquamarine', 'DeepSkyBlue', 'Aqua',
               'LightCyan', 'Turquoise', 'Aqua'] 
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    ax.spines['top'].set_color('w')
    ax.spines['bottom'].set_color('w')
    ax.spines['left'].set_color('w')
    ax.spines['right'].set_color('w')        
    ax.xaxis.label.set_color('w')
    ax.tick_params(axis='both', colors='w', length=8, width=4)
        
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(4)
    
    ct = 0
    #color = 0.2
    for j, name in enumerate(datasets):
    
        path = '/Users/lisalocey/Desktop/RareBio/data/macro'
        RADs = ff.get_SADs(path, name)
        
        if len(RADs) > 130:
            RADs = random.sample(RADs, 130)
        
        #r = lambda: random.randint(0,255)
        color = colors[j] # '#%02X%02X%02X' % (r(),r(),r())
        
        N = 0
        S = 0
        print name, len(RADs)
        
        for RAD in RADs:
            
            N = sum(RAD)
            S = len(RAD)
            
            RAD.sort()
            RAD.reverse()
        
            RAD = list([x for x in RAD if x != 0])
            #RAD = list(np.log(RAD))
        
            S = len(RAD)
        
            if S >= 20:
                ranks = range(1,len(RAD)+1)
	        for i, rank in enumerate(ranks):
        	    ranks[i] = rank/S
            
                plt.plot(ranks, RAD, color = str(color), lw = 0.75, alpha=0.8)
    
        plt.plot([0],[-1], c = color, lw = 2, label=labels[j]) #+'\nN='+str(int(N))+', S='+str(int(S)))
        #plt.xlim(1, 400)
        plt.yscale('log')
        #plt.ylim(0,12.0)
        leg = plt.legend(loc=1, prop={'size':18})
        leg.draw_frame(False)
            
        leg.get_frame().set_alpha(0.0)
        for text in leg.get_texts():
            text.set_color('w')
        
        #plt.ylabel('Abundance', color='w', fontsize=25)
        #plt.xlabel('Rank', color='w', fontsize=25)
        plt.tick_params(axis='both', which='major', labelsize=22)
    
        #color += 0.1    
        
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/macroRADs.png',
             dpi=600, transparent=True, bbox_inches = 'tight', pad_inches=0.03) 
    plt.close()
    
    return
    

def socioecoRADs():

    datasets = []
    for name in os.listdir('/Users/lisalocey/Desktop/XRADs/econ'): 
        datasets.append(name)                                 
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    labels = datasets
    
    ct = 0
    #color = 0.2
    for j, name in enumerate(datasets):
    
        path = '/Users/lisalocey/Desktop/XRADs/econ'
        RADs = ff.get_SADs(path, name)
    
        r = lambda: random.randint(0,255)
        color = '#%02X%02X%02X' % (r(),r(),r())
        color = 'b'
        print name, len(RADs)
        
        for RAD in RADs:
        
            RAD.sort()
            RAD.reverse()
        
            RAD = list([x for x in RAD if x != 0])
            RAD = list(np.log(RAD))
        
            S = len(RAD)
        
            if S >= 50:
                ranks = range(1,len(RAD)+1)
	        for i, rank in enumerate(ranks):
        	    ranks[i] = rank/S
            
                plt.plot(ranks, RAD, color = str(color), lw = 0.5, alpha=0.8)
    
        plt.plot([0],[-1], c = str(color), lw = 2, label=labels[j])
        plt.xlim(-0.02, 1.4)
        plt.ylim(0,22)
        leg = plt.legend(loc=1, prop={'size':3.6})
        leg.draw_frame(False)
        plt.ylabel('Abundance', color='w', fontsize=16)
        plt.xlabel('Rank', color='w', fontsize=16)
        plt.tick_params(axis='both', which='major', labelsize=22)
    
        #color += 0.1    
        
    plt.savefig('/Users/lisalocey/Desktop/econRADs.png',
             dpi=600, bbox_inches = 'tight', pad_inches=0.03) 
    plt.close()
    
    return
              
                            
                                          
#NvS()                                                                     
microRADs()
macroRADs()          
#socioecoRADs()
AgAndBCI()