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
rc("figure",figsize=(6, 6))

import pickle

import aftModel_v2 as model

wdir = '../data/'

# Get summary parameters for help in making T(A) curves
results = pickle.load(open('logOverviewResults.pickle','rb'))
obsArgsList, pFit, obsPred = results
obsArgs = obsArgsList[0]
obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,geothermalGradient = obsArgs

# summary data
df = pd.read_csv(wdir + 'aft_allData.csv')

# model results
log = pd.read_csv('logAll.csv')

# chisq(4) 0.05 and 0.5 confidence regions
#sel50 = (log['Q1'] < 3.36)
sel95 = (log['Q1'] < 9.49)

# set up a 3x3 layout
fig, ax = plt.subplots(nrows=3,ncols=3)

letter = ['A','B','C', 'D','E','F', 'G','H','I']

# step through models for each summary sample id (each age, getting older)
for i in range(9):
    row = int(i/3)
    col = np.mod(i,3)
    
    sel = (log['sample'] == i)
    
    log95 = log[sel & sel95]
    #log50 = log[sel & (log['Q1'] < 3.36)]
    # for debug, write to file
    # log95.to_csv('log95'+str(i)+'.csv',index=False)
    
    # arrays of all possible parameter values
    p95 = log95.iloc[:,5:12].values
    #p50 = log50.iloc[:,3:10].values
    
    # find envelopes of T(A) paths and plot
    if len(p95) > 0:    
        A95 = list()
        for j in range(len(p95)):
            A,T = model.AT_from_p_dA(p95[j],p_T,dT,Tsmooth)
            A95.append(A)
        A95 = np.asarray(A95)   
        A95min = A95.min(axis=0)
        A95max = A95.max(axis=0)
        #ax[row,col].fill_betweenx(T,A95min,A95max,color='pink')
        ax[row,col].fill_betweenx(T,A95min,A95max,color='lightgrey')
    '''
    if len(p50) > 0:
        A50 = list()
        for j in range(len(p50)):
            A,T = model.AT_from_p_dA(p50[j],p_T,dT,Tsmooth)
            A50.append(A)
        A50 = np.asarray(A50)   
        A50min = A50.min(axis=0)
        A50max = A50.max(axis=0)
        ax[row,col].fill_betweenx(T,A50min,A50max,color='steelblue')
    '''
    # plot best model curve for alpha = 0, 1 and 2
    
    logSel = log[sel & (log['alpha']==2)]
    p_dA = logSel[logSel['rss']==logSel['rss'].min()].values[0,5:12]
    A,T = model.AT_from_p_dA(p_dA,p_T,dT,Tsmooth)             
    ax[row,col].plot(A,T,c='blue',linestyle='--')
    
    logSel = log[sel & (log['alpha']==1)]
    p_dA = logSel[logSel['rss']==logSel['rss'].min()].values[0,5:12]
    A,T = model.AT_from_p_dA(p_dA,p_T,dT,Tsmooth)             
    ax[row,col].plot(A,T,c='red',linestyle=':')
    
    logSel = log[sel & (log['alpha']==0)]
    p_dA = logSel[logSel['rss']==logSel['rss'].min()].values[0,5:12]
    A,T = model.AT_from_p_dA(p_dA,p_T,dT,Tsmooth)             
    ax[row,col].plot(A,T,c='black')
    
    # Add label and age of observation
    ax[row,col].text(95,5,letter[i],weight='bold')
    txt = '{0:3.1f} \u00B1 {1:3.1f} Ma'.format(df['ageMean'].iloc[i],
                                               df['ageMeanErr'].iloc[i])
    ax[row,col].text(95,20,txt)
    
    # format look of plot
    if col == 0:
        ax[row,col].set(ylabel=r'Temperature ( $\degree$C)')
    if row == 2:
        ax[row,col].set(xlabel='Age (Ma)')
    ax[row,col].set(xlim=[0,100],ylim=[-10,150])
    ax[row,col].invert_xaxis()
    ax[row,col].invert_yaxis()
    ax[row,col].axvline(43,linestyle='--',color='steelblue')
    ax[row,col].axvline(53,linestyle='--',color='steelblue')
    
plt.tight_layout()

fig.savefig(wdir+'aftOverview_samples.jpg',dpi=300)
fig.savefig(wdir+'aftOverview_samples.pdf')

