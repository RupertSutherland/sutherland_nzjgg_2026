#!/usr/bin/env python
# python3
__author__ = "Rupert Sutherland"
"""
    load and plot summary data to model
"""
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import rc
rc("pdf", fonttype=42)
rc("font",size=6,family='Arial')
rc("figure",figsize=(5, 5))

wdir = '../data/'

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

fig.savefig(wdir+'aftSummary.jpg')
fig.savefig(wdir+'aftSummary.pdf')




