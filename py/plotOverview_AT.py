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

import pickle

wdir = '../data/'

results = pickle.load(open('logOverviewResults.pickle','rb'))
obsArgsList, pFit, obsPred = results

fig, ax = plt.subplots(3)

for i in range(0,3,1):
    
    ax[i].set(xlabel='Age (Ma)')
    ax[i].set(ylabel='Temperature ( $^\circ$C)')
    ax[i].set_xlim([0,100])
    ax[i].set_ylim([-10,140])
    ax[i].invert_xaxis()
    ax[i].invert_yaxis()
    ax[i].axvline(43,linestyle='--',color='steelblue')
    ax[i].axvline(53,linestyle='--',color='steelblue')
    
    for j in range(len(obsArgsList)):
        
        obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,geothermalGradient = obsArgsList[j]
        
        if alpha == i:
            
            A,T = model.AT_from_p_dA(pFit[j],p_T,dT,Tsmooth)
            
            ax[i].plot(A,T,color='grey')
            
            ax[i].text(95,30,r'$\alpha = $'+str(i))

ax[0].set_title('A',loc='left')
ax[1].set_title('B',loc='left')
ax[2].set_title('C',loc='left')

plt.tight_layout()

fig.savefig(wdir+'aftOverview_AT.jpg',dpi=300)
fig.savefig(wdir+'aftOverview_AT.pdf')

