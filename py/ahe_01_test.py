#!/usr/bin/env python
# python3
"""
    Produce T(t) for a target AFT age;
    check model aft age is correct;
    Compute sensitivity analysis for ahe age for range of parameters
    radius, U, Th
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc("pdf", fonttype=42)
rc("font",size=10,family='Arial')

from scipy.interpolate import interp1d

import aftModel_v3 as model

__author__ = "Rupert Sutherland"

obsArgs,model_output_bins = model.load_model_output_bins()
obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,geothermalGradient = obsArgs

aft_target = 35.6

p_dA = model.p_dA_from_aft(aft_target,obsArgs,model_output_bins)

t,T = model.tT_from_p_dA(p_dA,obsArgs)
#plt.plot(t,T)

# check AFT age
print('aft =', model.aft_from_tT(t, T, obsArgs)[0])

# -----------------------------
radius = 1e-3
U = 1e-5
Th = 1e-6

for radius in np.linspace(1e-4,1e-2,10):
#for U in np.linspace(1e-7,1e-4,10):
#for Th in np.linspace(1e-7,1e-4,10):

    ahe = model.ahe_from_tT(t,T,radius,U,Th)
    print('ahe =', ahe, 'radius,U,Th =',radius,U,Th)
