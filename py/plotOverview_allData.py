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

import pickle

wdir = '../data/'

#------------------------------------------
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
'''
ax[2].set(ylabel=r'Exhumation rate (m/Myr)')
ax[2].errorbar(df['ageMean'],df['exRate'],
               xerr=df['ageMeanErr'],yerr=df['exRateErr'], 
               fmt='.',markersize=2,color='black',ls='none',
               ecolor='grey',elinewidth=1)
'''
ax[2].set(ylabel=r'Age gradient (Myr/km)')
ax[2].errorbar(df['ageMean'],1000.*df['ageGrad'],
               xerr=df['ageMeanErr'],yerr=1000.*df['ageGradErr'], 
               fmt='.',markersize=2,color='black',ls='none',
               ecolor='grey',elinewidth=1)
ax[2].set(xlabel=r'Age (Ma)')

#------------------------------------------
results = pickle.load(open('logOverviewResults.pickle','rb'))
obsArgsList, pFit, obsPred = results

geotherm = obsArgsList[0][9]

alpha = np.asarray(obsArgsList)[:,4]

# Step through all results and plot alpha 0,1,2
col = ['black','steelblue','darkred']
style =['--','-',':']
for i in range(3):
    sel = [alpha == i]
    
    ax[0].set_title('A',loc='left')
    ax[0].plot(obsPred[sel][:,0],obsPred[sel][:,1],
               c=col[i],linestyle=style[i])
    ax[1].set_title('B',loc='left')
    ax[1].plot(obsPred[sel][:,0],obsPred[sel][:,2],
               c=col[i],linestyle=style[i])
    ax[2].set_title('C',loc='left')
    ax[2].plot(obsPred[sel][:,0],1000.*obsPred[sel][:,3]*geotherm,
               c=col[i],linestyle=style[i])

plt.tight_layout()

fig.savefig(wdir+'aftOverview_allData.jpg',dpi=300)
fig.savefig(wdir+'aftOverview_allData.pdf')

