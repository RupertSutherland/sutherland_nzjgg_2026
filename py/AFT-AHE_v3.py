#!/usr/bin/env python
# python3
"""
    Analyse cooling history model output for NVL...
    Load overview results and specific models.
    Specifically, the p_T temperature parameters (and others in obsArgs)
    For the 9 age bins defined (approx 5 Myr intervals), then find:
    aft_model_bins, the AFT model age constructed from parameters p_dA
    p_dA_bins, the parameters p_dA for each bin.
    
    
    
"""
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
from matplotlib import rc
rc("pdf", fonttype=42)
rc("font",size=10,family='Arial')

import aftModel_v2 as model

__author__ = "Rupert Sutherland"

obsArgs,model_output_bins = model.load_model_output_bins()
obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,geothermalGradient = obsArgs

aft_model_array = np.arange(20,71,1)
age_uniform_model_array = np.arange(0,81,1)

# for each age in array, compute T(A) and save
# Resample A onto uniform age_model_array
T_history_model_list = list()
for age in aft_model_array:
    p_dA = model.p_dA_from_aft(age,obsArgs,model_output_bins)
    A,T = model.AT_from_p_dA(p_dA,p_T,dT,Tsmooth)
    #plt.plot(-A,-T)
    interp_T_from_A = interp1d(A,T,fill_value='extrapolate')
    T_interp = interp_T_from_A(age_uniform_model_array)
    plt.plot(-age_uniform_model_array,-T_interp)
    
    '''
    # check age is correct
    aft_model = model.aft_from_p_dA(p_dA,obsArgs)
    print(age, age-aft_model[0])
    '''


