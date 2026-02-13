#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot overview log

@author: rupert
"""

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import rc
rc("pdf", fonttype=42)
rc("font",size=6,family='Arial')
rc("figure",figsize=(3, 6))

import aftModel_v2 as model

def strToArray(s):
    '''
    Converts a string output onto one line with str(array) i.e. of form
    '[[1 2][3 4]]' and returns np.array([1,2],[3,4])
    

    Parameters
    ----------
    s : str
        format '[[1 2][3 4]]'.
        Can be '[[1 2][3 4]]\n'

    Returns
    -------
    a : nd-array
        format np.array([1,2],[3,4]).

    '''
    from ast import literal_eval
    sSingleSpaced = ' '.join(s.split())
    sFormatted = sSingleSpaced.replace('[ ' , '[')
    sFormatted = sFormatted.replace(' ]' , ']')
    sFormatted = sFormatted.replace(' ' , ',')
    sFormatted = sFormatted.replace('][' , '],[')
    #print(sFormatted)
    return np.array(literal_eval(sFormatted))

#a = strToArray('[[1 2][3 4]]\n')

wdir = '../data/'


'''
Each block of logOverview.txt has format...
NEW_MODEL
obs
obsErr
Dpar
coolRate
alpha
p_T
dT
Tsmooth
Tdiff
geothermalGradient
qFit 0.0002880614479274299
pFit [22.94428836  5.2639607   0.08368587  0.08341437  0.0663763  10.80819383 21.06328781]
time(s): 117.42861342299966

obsArgs =  (obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,geothermalGradient)
obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,geothermalGradient = obsArgs
'''

obsArgsList = list()
p = list()
obsPred = list()

with open('logOverview.txt','r') as f:
    lines = f.readlines()

for i in range(len(lines)):
    if lines[i] == 'NEW_MODEL\n':
        obs = strToArray(lines[i+1])
        obsErr = strToArray(lines[i+2])
        Dpar = float(lines[i+3])
        coolRate = strToArray(lines[i+4])
        alpha = float(lines[i+5])
        p_T = strToArray(lines[i+6]) 
        dT = float(lines[i+7])
        Tsmooth = float(lines[i+8])
        Tdiff = float(lines[i+9])
        geothermalGradient = float(lines[i+10])
        pFit = strToArray(lines[i+12].replace('pFit ',''))
        
        obsArgs = (obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,
                        geothermalGradient)
        
        pred = model.aftGrad_from_p_dA(pFit,obsArgs)[0:4]
        
        obsArgsList.append(obsArgs)
        p.append(pFit)
        obsPred.append(pred)
        
p = np.asarray(p)
obsPred = np.asarray(obsPred)

fig, ax = plt.subplots(3)

for i in range(0,3,1):
    
    ax[i].set(xlabel='Age (Ma)')
    ax[i].set(ylabel='Temperature ( $^\circ$C)')
    ax[i].set_xlim([0,100])
    ax[i].set_ylim([-10,140])
    ax[i].invert_xaxis()
    ax[i].invert_yaxis()
    ax[i].axvline(43.5,linestyle='--',color='steelblue')
    ax[i].axvline(52.9,linestyle='--',color='steelblue')
    
    for j in range(len(obsArgsList)):
        
        obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,geothermalGradient = obsArgsList[j]
        
        if alpha == i:
            
            A,T = model.AT_from_p_dA(p[j],p_T,dT,Tsmooth)
            
            ax[i].plot(A,T,color='grey')
            
            ax[i].text(95,30,r'$\alpha = $'+str(i))

plt.tight_layout()

fig.savefig(wdir+'aftOverview_TA.jpg')
fig.savefig(wdir+'aftOverview_TA.pdf')

df = pd.read_csv(wdir + 'aft_allData.csv')

fig, ax = plt.subplots(3)

ax[0].set(ylabel=r'Mean track length ( $\mu$m)')
ax[0].errorbar(df['ageMean'],df['Lm_mean'],
               xerr=df['ageMeanErr'],yerr=df['Lm_err'], 
               fmt='.',markersize=2,color='black',ls='none',
               ecolor='grey',elinewidth=1)

ax[1].set(ylabel=r'Track length s.d. ( $\mu$m)')
ax[1].errorbar(df['ageMean'],df['sd_mean'],
               xerr=df['ageMeanErr'],yerr=df['sd_err'], 
               fmt='.',markersize=2,color='black',ls='none',
               ecolor='grey',elinewidth=1)

ax[2].set(ylabel=r'Exhumation rate (m/Myr)')
ax[2].errorbar(df['ageMean'],df['exRate'],
               xerr=df['ageMeanErr'],yerr=df['exRateErr'], 
               fmt='.',markersize=2,color='black',ls='none',
               ecolor='grey',elinewidth=1)

plt.tight_layout()

fig.savefig(wdir+'aftOverview_allData.jpg')
fig.savefig(wdir+'aftOverview_allData.pdf')


