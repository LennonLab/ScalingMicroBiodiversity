from __future__ import division
import  matplotlib.pyplot as plt

import pandas as pd
import linecache
import numpy as np
import os
import sys

#import statsmodels.stats.api as sms
#import statsmodels.api as sm
import statsmodels.formula.api as smf
#from statsmodels.sandbox.regression.predstd import wls_prediction_std
from statsmodels.stats.outliers_influence import summary_table

mydir = os.path.expanduser("~/GitHub/MicrobialScaling/")
mydir2 = os.path.expanduser("~/")


def Fig1(condition, ones, sampling):

    tail = str()
    if ones is False:
        tail = '-SADMetricData_NoMicrobe1s.txt'
    elif ones is True:
        tail = '-SADMetricData.txt'

    datasets = []
    GoodNames = []
    emp = str()

    if condition == 'open': emp = 'EMPopen'
    elif condition == 'closed': emp = 'EMPclosed'

    #GoodNames = [emp, 'HMP', 'BBS', 'CBC', 'MCDB', 'GENTRY', 'FIA']
    #GoodNames = [emp, 'BBS', 'CBC', 'MCDB', 'GENTRY', 'FIA'] # all microbe data is emp
    #GoodNames = ['HMP', 'BBS', 'CBC', 'MCDB', 'GENTRY', 'FIA'] # all microbe data is HMP
    GoodNames = [emp, 'HMP', 'BIGN', 'TARA', 'BOVINE', 'HUMAN', 'LAUB', 'SED', 'CHU', 'CHINA', 'CATLIN', 'FUNGI', 'HYDRO', 'BBS', 'CBC', 'MCDB', 'GENTRY', 'FIA'] # all microbe data is MGRAST

    for name in os.listdir(mydir +'data/micro'):
        if name in GoodNames: pass
        else: continue

        path = mydir+'data/micro/'+name+'/'+name+tail
        num_lines = sum(1 for line in open(path))
        datasets.append([name, 'micro', num_lines])
        print name, num_lines


    for name in os.listdir(mydir +'data/macro'):
        if name in GoodNames: pass
        else: continue

        path = mydir+'data/macro/'+name+'/'+name+'-SADMetricData.txt'
        num_lines = sum(1 for line in open(path))
        datasets.append([name, 'macro', num_lines])
        print name, num_lines


    metrics = ['Rarity, '+r'$log_{10}$',
            'Dominance, '+r'$log_{10}$',
            'Evenness, ' +r'$log_{10}$',
            'Richness, ' +r'$log_{10}$',] #+r'$(S)^{2}$']

    fig = plt.figure()
    for index, i in enumerate(metrics):

        metric = i
        fig.add_subplot(2, 2, index+1)
        fs = 12 # font size used across figures

        MicIntList, MicCoefList, MacIntList, MacCoefList, R2List, metlist = [[], [], [], [], [], []]
        Nlist, Slist, Evarlist, ESimplist, klist, radDATA, BPlist, NmaxList, rareSkews, KindList, StdList = [[], [], [], [], [], [], [], [], [], [], []]
        #name, kind, N, S, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]

        its = 10

        for n in range(its):

            #name, kind, N, S, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
            Nlist, Slist, Evarlist, ESimplist, klist, radDATA, BPlist, NmaxList, rareSkews, KindList, StdList = [[], [], [], [], [], [], [], [], [], [], []]

            numMac = 0
            numMic = 0
            radDATA = []

            for dataset in datasets:

                name, kind, numlines = dataset
                lines = []
                small = ['BIGN', 'BOVINE', 'CHU', 'LAUB', 'SED']
                big = ['HUMAN', 'CHINA', 'CATLIN', 'FUNGI', 'HYDRO']

                if kind == 'macro':
                    lines = np.random.choice(range(1, numlines+1), 100, replace=True)

                elif name in small:
                    lines = np.random.choice(range(1, numlines+1), 20, replace=True)

                elif name in big:
                    lines = np.random.choice(range(1, numlines+1), 50, replace=True)

                elif name == 'TARA':
                    lines = np.random.choice(range(1, numlines+1), 50, replace=True)
                else:
                    lines = np.random.choice(range(1, numlines+1), 50, replace=True)

                if kind == 'micro':
                    path = mydir+'data/'+kind+'/'+name+'/'+name+tail

                else:
                    path = mydir+'data/macro/'+name+'/'+name+'-SADMetricData.txt'

                for line in lines:
                    data = linecache.getline(path, line)
                    radDATA.append(data)

            for data in radDATA:

                data = data.split()
                name, kind, N, S, Var, Evar, ESimp, EQ, O, ENee, EPielou, EHeip, BP, SimpDom, Nmax, McN, skew, logskew, chao1, ace, jknife1, jknife2, margalef, menhinick, preston_a, preston_S = data

                N = float(N)
                S = float(S)

                Nlist.append(float(np.log10(N)))
                Slist.append(float(np.log10(S)))

                ESimplist.append(float(np.log10(float(ESimp))))
                KindList.append(kind)

                BPlist.append(float(BP))
                NmaxList.append(float(np.log10(float(Nmax))))

                # log-modulo transformation of skewnness
                lms = np.log10(np.abs(float(skew)) + 1)
                if skew < 0: lms = lms * -1
                rareSkews.append(float(lms))

                if kind == 'micro':
                    numMic += 1
                    klist.append('b')
                if kind == 'macro':
                    klist.append('r')
                    numMac += 1

            if index == 0: metlist = list(rareSkews)
            elif index == 1: metlist = list(NmaxList)
            elif index == 2: metlist = list(ESimplist)
            elif index == 3: metlist = list(Slist)

            # Multiple regression
            d = pd.DataFrame({'N': list(Nlist)})
            d['y'] = list(metlist)
            d['Kind'] = list(KindList)

            f = smf.ols('y ~ N * Kind', d).fit()
            #f = smf.rlm('y ~ N * Kind', d).fit()
            #r2 = smf.wls('y ~ N * Kind', d, weights= f.weights).fit().rsquared
            r2 = f.rsquared

            MacIntList.append(f.params[0])
            MacCoefList.append(f.params[2])

            if f.pvalues[1] < 0.05:
                MicIntList.append(f.params[1] + f.params[0])
            else:
                MicIntList.append(f.params[0])

            if f.pvalues[3] < 0.05:
                MicCoefList.append(f.params[3] + f.params[2])
            else:
                MicCoefList.append(f.params[2])

            R2List.append(r2)


        MacPIx, MacFitted, MicPIx, MicFitted = [[],[],[],[]]
        macCiH, macCiL, micCiH, micCiL = [[],[],[],[]]

        MacListX = []
        MacListY = []
        MicListX = []
        MicListY = []

        for j, k in enumerate(KindList):
            if k == 'micro':
                MicListX.append(Nlist[j])
                MicListY.append(metlist[j])

            elif k == 'macro':
                MacListX.append(Nlist[j])
                MacListY.append(metlist[j])

        print metric

        ols = smf.ols('y ~ N * Kind', d).fit()

        #rlm = smf.rlm('y ~ N * Kind', d).fit()
        #wls = smf.wls('y ~ N * Kind', d, weights= rlm.weights).fit()
        #r2 = wls.rsquared
        r2 = ols.rsquared

        st, data, ss2 = summary_table(ols, alpha=0.05)
        # ss2: Obs, Dep Var Population, Predicted Value, Std Error Mean Predict,
        # Mean ci 95% low, Mean ci 95% upp, Predict ci 95% low, Predict ci 95% upp,
        # Residual, Std Error Residual, Student Residual, Cook's D

        #fittedvalues = data[:,2]
        #predict_mean_se = data[:,3]
        predict_mean_ci_low, predict_mean_ci_upp = data[:,4:6].T
        predict_ci_low, predict_ci_upp = data[:,6:8].T

        for j, kval in enumerate(KindList):
            if kval == 'macro':

                macCiH.append(predict_mean_ci_upp[j])
                macCiL.append(predict_mean_ci_low[j])
                MacPIx.append(Nlist[j])
                MacFitted.append(ols.fittedvalues[j])

            elif kval == 'micro':

                micCiH.append(predict_mean_ci_upp[j])
                micCiL.append(predict_mean_ci_low[j])
                MicPIx.append(Nlist[j])
                MicFitted.append(ols.fittedvalues[j])

        MicPIx, MicFitted, micCiH, micCiL = zip(*sorted(zip(MicPIx, MicFitted, micCiH, micCiL)))
        MacPIx, MacFitted, macCiH, macCiL = zip(*sorted(zip(MacPIx, MacFitted, macCiH, macCiL)))

        num = min(len(MacListX), len(MicListX))
        micnums = np.random.choice(range(0, len(MicListX)), num, replace=False)
        macnums = np.random.choice(range(0, len(MacListX)), num, replace=False)


        for i, ind in enumerate(micnums):
            plt.scatter(MacListX[macnums[i]], MacListY[macnums[i]], color = 'LightCoral', alpha= 1 , s = 8, linewidths=0.5, edgecolor='Crimson')
            plt.scatter(MicListX[ind], MicListY[ind], color = 'SkyBlue', alpha= 1 , s = 8, linewidths=0.5, edgecolor='Steelblue')

        plt.fill_between(MacPIx, macCiL, macCiH, color='LightCoral', lw=0.0, alpha=0.9)
        plt.plot(MacPIx, MacFitted,  color='r', ls='--', lw=0.5, alpha=0.9)
        plt.fill_between(MicPIx, micCiL, micCiH, color='b', lw=0.0, alpha=0.3)
        plt.plot(MicPIx, MicFitted,  color='b', ls='--', lw=0.5, alpha=0.9)

        MicInt = round(np.mean(MicIntList), 2)
        MicCoef = round(np.mean(MicCoefList), 2)
        MacInt = round(np.mean(MacIntList), 2)
        MacCoef = round(np.mean(MacCoefList), 2)
        R2 = round(np.mean(R2List), 2)

        if index == 0:
            plt.ylim(-0.1, 2.5)
            plt.xlim(0, 8.2)

            plt.text(0.35, 2.1, r'$micro$'+ ' = '+str(round(10**MicInt,2))+'*'+r'$N$'+'$^{'+str(round(MicCoef,2))+'}$', fontsize=fs, color='Steelblue')
            plt.text(0.35, 1.8, r'$macro$'+ ' = '+str(round(10**MacInt,2))+'*'+r'$N$'+'$^{'+str(round(MacCoef,2))+'}$', fontsize=fs, color='Crimson')
            plt.text(0.35, 1.4,  r'$R^2$' + '=' +str(R2), fontsize=fs-1, color='k')

            plt.scatter([0],[-1], color = 'SkyBlue', alpha = 1, s=15, linewidths=0.9, edgecolor='Steelblue', label= 'microbes (n='+str(len(MicListY))+')')
            plt.scatter([0],[-1], color = 'LightCoral',alpha= 1, s=15, linewidths=0.9, edgecolor='Crimson', label= 'macrobes (n='+str(len(MacListY))+')')
            plt.legend(bbox_to_anchor=(-0.04, 1.05, 2.48, .2), loc=10, ncol=2, mode="expand",prop={'size':fs})

        elif index == 1:

            plt.plot([0,8.2],[0,8.2], ls = '--', lw=1, c='0.7')
            plt.ylim(0, 8)
            plt.xlim(0, 8.2)

            plt.text(0.35, 6.7, r'$micro$'+ ' = '+str(round(10**MicInt,2))+'*'+r'$N$'+'$^{'+str(round(MicCoef,2))+'}$', fontsize=fs, color='Steelblue')
            plt.text(0.35, 5.7, r'$macro$'+ ' = '+str(round(10**MacInt,2))+'*'+r'$N$'+'$^{'+str(round(MacCoef,2))+'}$', fontsize=fs, color='Crimson')
            plt.text(0.35, 4.7,  r'$R^2$' + '=' +str(R2), fontsize=fs-1, color='k')

        elif index == 2:
            plt.ylim(-3.5, 0.0)
            plt.xlim(0, 8.2)

            plt.text(0.35, -2.9, r'$micro$'+ ' = '+str(round(10**MicInt,2))+'*'+r'$N$'+'$^{'+str(round(MicCoef,2))+'}$', fontsize=fs, color='Steelblue')
            plt.text(0.35, -3.3, r'$macro$'+ ' = '+str(round(10**MacInt,2))+'*'+r'$N$'+'$^{'+str(round(MacCoef,2))+'}$', fontsize=fs, color='Crimson')
            plt.text(0.35, -2.5,  r'$R^2$' + '=' +str(R2), fontsize=fs-1, color='k')

        elif index == 3:
            plt.ylim(0.9, 5.0)
            plt.xlim(0, 8.2)

            plt.text(0.35, 4.5, r'$micro$'+ ' = '+str(round(2**MicInt,2))+'*'+r'$N$'+'$^{'+str(round(MicCoef,2))+'}$', fontsize=fs, color='Steelblue')
            plt.text(0.35, 4.0, r'$macro$'+ ' = '+str(round(2**MacInt,2))+'*'+r'$N$'+'$^{'+str(round(MacCoef,2))+'}$', fontsize=fs, color='Crimson')
            plt.text(0.35, 3.5,  r'$R^2$' + '=' +str(R2), fontsize=fs-1, color='k')
            print condition, ones, ': S =', '%.3e' % (10**(MicInt + MicCoef*(30.0)))
            #print R2


        plt.xlabel('$log$'+r'$_{10}$'+'($N$)', fontsize=fs)
        plt.ylabel(metric, fontsize=fs)
        plt.tick_params(axis='both', which='major', labelsize=fs-3)

    plt.subplots_adjust(wspace=0.4, hspace=0.4)

    if ones == False:
        plt.savefig(mydir+'/figs/Fig1/Locey_Lennon_2015_Fig1-'+condition+'_NoSingletons_'+str(sampling)+'.pdf', dpi=300, bbox_inches = "tight")
    if ones == True:
        plt.savefig(mydir+'/figs/Fig1/Locey_Lennon_2015_Fig1-'+condition+'_'+str(sampling)+'.pdf', dpi=300, bbox_inches = "tight")

    #plt.show()
    plt.close()

    return



EMPcondition = ['closed']
Singletons = [False]
#Samplings = [50, 100, 200, 300, 400, 500, 1000, 5000, 15000]
Samplings = [100]

for condition in EMPcondition:
    for ones in Singletons:
        for sampling in Samplings:
            Fig1(condition, ones, sampling)
