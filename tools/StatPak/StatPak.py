# -*- coding: utf-8 -*-
from __future__ import division
import  matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
import random
import scipy as sc
import os
import sys

import pandas
from pandas.tools import plotting
from scipy import stats
import statsmodels
from statsmodels.formula.api import ols
from numpy.random import randn


""" http://statsmodels.sourceforge.net/devel/stats.html#residual-diagnostics-and-specification-tests """


def NormStats(resids):
    
    DW = statsmodels.stats.stattools.durbin_watson(resids, axis=0) # Calculate the Durbin-Watson statistic for normality

    JB = statsmodels.stats.stattools.jarque_bera(resids, axis=0) # Calculate the Jarge-Bera test 

    Omni = statsmodels.stats.stattools.omni_normtest(resids, axis=0) # Calculate the Omnibus test for normal skewnness and kurtosis

    NormAd = statsmodels.stats.diagnostic.normal_ad(x, axis=0) # Anderson-Darling test for normal distribution unknown mean and variance

    KSnorm = statsmodels.stats.diagnostic.kstest_normal(x, pvalmethod='approx') # Lillifors test for normality, Kolmogorov Smirnov test with estimated mean and variance

    Lfor = statsmodels.stats.diagnostic.lillifors(x, pvalmethod='approx') # Lillifors test for normality, Kolmogorov Smirnov test with estimated mean and variance

    return

def AutoCorrStats(x, results, lags=None, nlags=None, store=False, boxpierc=False):

    Lj = statsmodels.stats.diagnostic.acorr_ljungbox(x, lags=None, boxpierce=False) # Calculate the Ljung-Box test for no autocorrelation

    BG = statsmodels.stats.diagnostic.acorr_breush_godfrey(results, nlags=None, store=False) # Calculate the Breush Godfrey Lagrange Multiplier tests for residual autocorrelation

    return 


    

def ResidStat(resid, exog_het):

    HB = statsmodels.stats.diagnostic.het_breushpagan(resid, exog_het) # The tests the hypothesis that the residual variance does not depend on the variables in x in the form

    return

def HetScedStats(resid, exog):

    HW = statsmodels.stats.diagnostic.het_white(resid, exog, retres=False) # White’s Lagrange Multiplier Test for Heteroscedasticity

    HA = statsmodels.stats.diagnostic.het_arch(resid, maxlag=None, autolag=None, store=False, regresults=False, ddof=0) # Engle’s Test for Autoregressive Conditional Heteroscedasticity (ARCH)

    return


def Linearity(res, resid, exog, olsresidual, olsresults):

    LH = statsmodels.stats.diagnostic.linear_harvey_collier(res) # Harvey Collier test for linearity. The Null hypothesis is that the regression is correctly modeled as linear.

    LR = statsmodels.stats.diagnostic.linear_rainbow(res, frac=0.5) # Rainbow test for linearity, The Null hypothesis is that the regression is correctly modelled as linear. The alternative for which the power might be large are convex, check.
 
    Llm = statsmodels.stats.diagnostic.linear_lm(resid, exog, func=None) # Lagrange multiplier test for linearity against functional alternative

    Bcusum = statsmodels.stats.diagnostic.breaks_cusumolsresid(olsresidual, ddof=0) # cusum test for parameter stability based on ols residuals

    BH = statsmodels.stats.diagnostic.breaks_hansen(olsresults) # test for model stability, breaks in parameters for ols, Hansen 1992

    Rols = statsmodels.stats.diagnostic.recursive_olsresiduals(olsresults, skip=None, lamda=0.0, alpha=0.95) # calculate recursive ols with residuals and cusum test statistic


def Outliers(results):
    
    #class statsmodels.stats.outliers_influence.OLSInfluence(results)
    OutInf = statsmodels.stats.outliers_influence.OLSInfluence(results)

    return
    
""" Incomplete 
# HGQ = statsmodels.stats.diagnostic.HetGoldfeldQuandt # function is not complete in statsmodels documentation
# ComCox = class statsmodels.stats.diagnostic.CompareCox # Cox Test for non-nested models
"""

