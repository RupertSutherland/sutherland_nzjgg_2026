#!/usr/bin/env python
# python3
"""
    Analyse aft tracks
"""
__author__ = "Rupert Sutherland"

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib.pyplot as plt
from matplotlib import rc
rc("pdf", fonttype=42)
rc("font",size=7,family='Arial')
rc("figure",figsize=(5, 3))

def binY(x,y,xbins=10,xrange=None,statistic='mean'):
    """
    Finds statistical value of y values in each x bin, given (x,y) data. 
    See scipy.stats.binned_statistic() for options.
    
    Parameters
    ----------
    x : array
        x values.
    y : array
        y values.
    statistic : string or callable, optional
        See documentation for scipy.stats.binned_statistic(). Default is mean.
        'mean','std','median','count','sum','min','max', or user-defined func.
    xbins : int or sequence of scalars, optional
        If xbins is an integer, it is the number of equal bins within xrange.
        If xbins is an array, then it is the location of xbin edges, similar
        to definitions used by np.histogram. Default is 10 bins.
        All but the last (righthand-most) bin is half-open. In other words, if 
        bins is [1, 2, 3, 4], then the first bin is [1, 2) (including 1, but 
        excluding 2) and the second [2, 3). The last bin, however, is [3, 4], 
        which includes 4.    
        
    xrange : (float, float) or [(float, float)], optional
        The lower and upper range of the bins. If not provided, range is 
        simply (x.min(), x.max()). Values outside the range are ignored.
    
    Returns
    -------
    y_stat : array
        The y statistic (e.g. mean) in each bin.       
    """
    y_stat, xbin_edges, binnumber = stats.binned_statistic(x, y, 
                                 statistic=statistic, bins=xbins, range=xrange)
    return y_stat


wdir = '../data/'

file1 = wdir + 'aft_tam.csv'
file2 = wdir + 'ahe_tam.csv'

aft = pd.read_csv(file1)

aft = aft[aft['Lm'] > 0.]
aft = aft[aft['lat'] > -76]

ageMin = 18.0
ageMax = 85.0
ageInc = 10.0
ageBins = np.arange(ageMin,ageMax+ageInc,ageInc)

age_mean = binY(aft['age'],aft['age'],xbins=ageBins,statistic='mean')
n = binY(aft['age'],aft['age'],xbins=ageBins,statistic='count')
age_std  = binY(aft['age'],aft['age'],xbins=ageBins,statistic='std')
age_serr = age_std / np.sqrt(n - 1)

Lm_mean = binY(aft['age'],aft['Lm'],xbins=ageBins,statistic='mean')
Lm_std  = binY(aft['age'],aft['Lm'],xbins=ageBins,statistic='std')
Lm_serr = Lm_std / np.sqrt(n - 1)

sd_mean = binY(aft['age'],aft['L_sd'],xbins=ageBins,statistic='mean')
sd_std  = binY(aft['age'],aft['L_sd'],xbins=ageBins,statistic='std')
sd_serr = sd_std / np.sqrt(n - 1)

df = pd.DataFrame(
        data = np.column_stack((age_mean,age_serr,n,Lm_mean,Lm_serr,
                                sd_mean,sd_serr)),
        columns = ['age','age_err','n','Lm_mean','Lm_err','sd_mean','sd_err'])
                 
df.to_csv(wdir + 'aft_nvl_tracks.csv', index=False)                 
                  
fig, ax = plt.subplots(2)

ax[0].errorbar(age_mean,Lm_mean,xerr=age_serr,yerr=Lm_serr,ls='none')
ax[0].set(ylabel=r'Mean track length ( $\mu$m)')

ax[1].errorbar(age_mean,sd_mean,xerr=age_serr,yerr=sd_serr,ls='none')
ax[1].set(xlabel='AFT Age (Ma)',ylabel=r'Track length s.d. ( $\mu$m)')

plt.tight_layout()

fig.savefig(wdir+'aft_tracks_nvl.jpg')
fig.savefig(wdir+'aft_tracks_nvl.pdf')

print(n.sum(),n)
