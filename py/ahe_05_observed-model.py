#!/usr/bin/env python
# python3
"""
    Produce T(t) for a target AFT age;
    check model aft age is correct;
    Compute sensitivity analysis for ahe age for range of parameters
    radius, U, Th
    Plot results for radius
"""
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import rc
rc("pdf", fonttype=42)
rc("font",size=10,family='Arial')

from scipy.interpolate import interp1d
from scipy.stats import norm
from scipy.stats import normaltest

import aftModel_v3 as model

__author__ = "Rupert Sutherland"

def cdf(xSorted):
    '''
    Takes a sorted (increasing) 1D array and returns corresponding 
    cumulative frequency values for each point (bin-centred).

    Parameters
    ----------
    xSorted : 1d array
        Input data.

    Returns
    -------
    1d array
        Output cdf values for each value of xSorted.

    '''
    n = len(xSorted)
    df = 100./(2*n)
    return np.linspace(df,100-df,n)

'''
'''

obsArgs,model_output_bins = model.load_model_output_bins()
obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,geothermalGradient = obsArgs

# load observations
data = pd.read_csv('../data/aft-ahe_data.csv')
aft = data['age_AFT'].values
ahe = data['age_AHE'].values
n = len(aft)

# Setup figure
plt.figure(figsize=(4,3))
xmin = -30
xmax = 200
plt.xlim(xmin,xmax)
plt.xlabel('AHE observed - expected (Myr)') 
plt.ylabel(r'Cumulative frequency (%)')
plt.minorticks_on()
plt.plot([0,0],[0,100],linestyle='-',c='grey')

# find AHE model age for each AHE observation, using T(t) that fits AFT age.
ahe_model = np.zeros_like(ahe)  
for i in range(n):
    p_dA = model.p_dA_from_aft(aft[i],obsArgs,model_output_bins)   
    t,T = model.tT_from_p_dA(p_dA,obsArgs)
    ahe_model[i] = model.ahe_from_tT(t,T,radius=500e-6,U=50e-6,Th=10e-6)              
ageDiff = ahe - ahe_model
ageDiffSorted = np.sort(ageDiff)
cdf_ageDiffSorted = cdf(ageDiffSorted)

plt.scatter(ageDiffSorted,cdf_ageDiffSorted,
            marker='o',c='steelblue',label='R = 500 $\mu$m')

# find AHE model age for each AHE observation, using T(t) that fits AFT age.
ahe_model = np.zeros_like(ahe)  
for i in range(n):
    p_dA = model.p_dA_from_aft(aft[i],obsArgs,model_output_bins)   
    t,T = model.tT_from_p_dA(p_dA,obsArgs)
    ahe_model[i] = model.ahe_from_tT(t,T,radius=50e-6,U=50e-6,Th=10e-6)              
ageDiff2 = ahe - ahe_model
ageDiffSorted2 = np.sort(ageDiff2)
cdf_ageDiffSorted2 = cdf(ageDiffSorted2)

plt.scatter(ageDiffSorted2,cdf_ageDiffSorted2,
            marker='+',c='darkred',label='R = 50 $\mu$m')

plt.legend()
plt.tight_layout()
plt.savefig('ahe_obs-model.jpg',dpi=300)
plt.savefig('ahe_obs-model.pdf')

