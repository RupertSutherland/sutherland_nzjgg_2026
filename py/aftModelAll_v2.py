#!/usr/bin/env python
# python3
__author__ = "Rupert Sutherland"
"""
    Takes 3.4 hr to run on linux vm with f90 compiled and dT = 1.
"""
import numpy as np
import pandas as pd

import time
from shutil import copy

import aftModel_v2 as model

Dpar = 1.7

tStart = time.process_time()

wdir = '../data/'

df = pd.read_csv(wdir + 'aft_allData.csv')

df['coolRate'] = abs(df['exRate'] * model.geothermalGradient)
df['coolRateErr'] = abs(df['exRateErr'] * model.geothermalGradient)
coolRateList = ['ageMean','coolRate','coolRateErr']
coolRate = df[coolRateList].values

with open('logOverview.txt','w') as logOverview:
    logOverview.write('obs\nobsErr\nDpar\ncoolRate\nalpha\np_T\ndT\nTsmooth\n'+
                      'Tdiff\ngeothermalGradient\n\n')
    
# Overall Fit is Q = Q1 + Q2 * (alpha / N)
# alpha is weight given to rss of N global cooling rate values Q2
# Q1 is rss of fit to [age,Lmean,Lsd,ageGrad_m] value(s)
for alpha in range(10):
                      
    for i in range(len(df)):
    #for i in range(2):
        
        iStr = '{:02d}'.format(i)
        
        #-------------------------------------------------------------------------
        # Array of observations [obs0,obs1...] where obs = [age,Lmean,Lsd,ageGrad_m]
        # Predicted ageGrad is ageGrad_T in Myr/C (negative sign)
        # To convert to Myr/m (spatial, positive)
        # dA/dz = dA/dT * dT/dz
        # But z positive up, T positive down; assume geothermalGradient = -0.03 C/m
        # ageGrad_m = - ageGrad_T * geothermalGradient
        # exhumation at 1 mm/yr = 1000 m/Myr gives dA/dz = 0.001 Myr/m
        
        #[age,Lmean,Lsd,ageGrad_m]
        obs =    np.array([df['ageMean'].iloc[i], 
                           df['Lm_mean'].iloc[i], 
                           df['sd_mean'].iloc[i], 
                           df['ageGrad'].iloc[i] 
                           ]) 
        obsErr = np.array([df['ageMeanErr'].iloc[i], 
                           df['Lm_err'].iloc[i], 
                           df['sd_err'].iloc[i], 
                           df['ageGradErr'].iloc[i] 
                           ]) 
           
        obsArgs = (obs,obsErr,Dpar,coolRate,alpha,
                   model.p_T, model.dT, model.Tsmooth, model.Tdiff,
                   model.geothermalGradient)
        # obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,geothermalGradient = obsArgs
        
        with open('logOverview.txt','a') as logOverview:
            logOverview.write('NEW_MODEL\n')
            for v in obsArgs:
                logOverview.write(str(v).replace('\n','')+'\n')
        
        # Guess first parameters, uniform cooling
        # Ages at reference temperatures p_T
        p_A = obs[0] * model.p_T / (model.Tclosure - model.p_T[0])
        # Age differences - our choice of adjustable parameter
        p_dA = np.diff(p_A)
        # Bounding values of dA to control optimization
        pBounds = np.column_stack([0.01*p_dA,100.*p_dA])
        
        model.logCreate(p_dA)  
        
        #qFit,pFit = (0,0)
        qFit,pFit = model.fitParameters(p_dA,obsArgs,pBounds)
           
        copy('logfile.txt', 'logfile' + str(alpha) + iStr + '.txt')
    
        #model.plotLogfile()
        
        with open('logOverview.txt','a') as logOverview:
            logOverview.write('qFit '+str(qFit)+
                              '\npFit '+str(pFit).replace('\n','')+'\n')
            logOverview.write('time(s): '+str(time.process_time() - tStart)+'\n\n')
        
        print('time',time.process_time() - tStart)
