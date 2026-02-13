#!/usr/bin/env python
# python3
"""
    Load and plot summary data
"""
__author__ = "Rupert Sutherland"

import numpy as np
from scipy.interpolate import interp1d

import pandas as pd


wdir = '../data/'

tracks = pd.read_csv(wdir + 'aft_nvl_tracks.csv')

# this file is finer time-sampled and includes exhumation rate
dat = pd.read_csv(wdir + 'aftGrad_nvl_binSummary.csv')

dat = dat[dat['ageMean']<85]

def resample(c):
    # resample tracks data onto ageMean of dat, which includes exhumation rate
    # c is the column name
    interpFunc = interp1d(tracks['age'],tracks[c],fill_value='extrapolate')
    newVals = interpFunc(dat['ageMean'])
    dat[c] = newVals
    
for c in ['Lm_mean','Lm_err','sd_mean','sd_err']:
    resample(c)

dat.to_csv(wdir + 'aft_allData.csv', index=False)
