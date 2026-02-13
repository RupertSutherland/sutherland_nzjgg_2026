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
import matplotlib.pyplot as plt
from matplotlib import rc
rc("pdf", fonttype=42)
rc("font",size=10,family='Arial')

from scipy.interpolate import interp1d

import aftModel_v3 as model

__author__ = "Rupert Sutherland"

obsArgs,model_output_bins = model.load_model_output_bins()
obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,geothermalGradient = obsArgs

U = 3e-5
Th = 1e-5
logRadius = np.linspace(-4.5,-2.5,50)
radius = 10**logRadius
ahe = np.zeros_like(radius)
    
plt.figure(figsize=(4,3))
plt.xscale('log')
plt.xlabel('Radius (mm)')
plt.ylabel('AFT-AHE (Ma)')

aftRange = np.arange(30,71,5)
n = len(aftRange)
colors = plt.cm.jet(np.linspace(1,0,n))
'''
'''   
aheList = list()
for aft in aftRange:
    print(aft)
    p_dA = model.p_dA_from_aft(aft,obsArgs,model_output_bins)   
    t,T = model.tT_from_p_dA(p_dA,obsArgs)
    #plt.plot(t,T)    
    # check AFT age
    aft_model = model.aft_from_tT(t, T, obsArgs)[0]
    print('aft =', aft_model) 
    
    ahe = np.zeros_like(radius)
    for j in range(len(ahe)): 
        
        ahe[j] = model.ahe_from_tT(t,T,radius[j],U,Th)  
        
    aheList.append(ahe)
    
for i in range(n):
    
    plt.plot(radius*1e3, aftRange[i]-aheList[i],color=colors[i],)

plt.minorticks_on()
plt.tight_layout()
plt.savefig('ahe_radius_sensitivity.jpg',dpi=300)
plt.savefig('ahe_radius_sensitivity.pdf')

