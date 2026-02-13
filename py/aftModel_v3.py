#!/usr/bin/env python
# python3
"""
    Model cooling histories to fit AFT data.
    Improve definition of p95-p50 confidence interval
    
    p_dA:
        Adjustable parameters. Array age differences [... dAi ...] 0<=i<n
        Constrain 0 < dAi < some geologically reasonable value, e.g. 10 mm/yr
    p_A:
        Derived. Adjustable parameters expressed as age. 
        The first age is 0.0 Ma, so 
        len(p_A) = len(p_dA) + 1 = n + 1
    p_T:
        Array of temperatures for each p_A.

    Linear fit through points, smoothed with sliding Gaussian window with
    width sigma = Tsmooth.   
    
"""
__author__ = "Rupert Sutherland"

import numpy as np
import pandas as pd

import time
import pickle
import itertools 
from scipy.interpolate import interp1d
from scipy.ndimage.filters import gaussian_filter1d
from scipy.optimize import minimize

import AFTannealingLib
import helium_diffusion_models

import matplotlib.pyplot as plt

#-------------------------------------------------------------------------
# GLOBAL PARAMETERS USED IN SOME FUNCTIONS
#reference temperatures, T increment, and smoothing factor
p_T = np.array([-10,0,30,60,90,110,130,160])
dT = 1.      # Temperature mesh spacing
Tsmooth = 5. # for smoothing T(t) paths
Tdiff = 5.  # finite value for computing dPred/dTdiff gradients
geothermalGradient = -0.03 # degC/m note dT/dz has z positive upwards
Tclosure = 100. # Guess at AFT closure temperature

def A_from_p_A(p_A,T,dT,Tsmooth):
    '''
    Compute age array from age parameters.

    Parameters
    ----------
    p_A : array-like
        Age parameters (not dA).
    T : nd-array
        Temperature array.
    dT : float
        Temperature increment value.
    Tsmooth : float
        Temperature smoothing window width 
        (for sigma of Gaus smoothing window).

    Returns
    -------
    Asmooth : nd-array
        Values of age at each temperature in T.

    '''
    A_interpFunc = interp1d(p_T, p_A, kind='linear', fill_value="extrapolate")
    nSmooth = int(Tsmooth/dT)
    # pad T with same gradient before and after by plenty so ends not smoothed
    Tpad = np.arange(dT,dT*nSmooth*5,dT)
    n = len(Tpad)
    Tpadded = np.concatenate([T[0] - np.flip(Tpad), T, T[-1] + Tpad])
    A = A_interpFunc(Tpadded)
    AsmoothPadded = gaussian_filter1d(A, sigma=nSmooth, mode='nearest') 
    Asmooth = AsmoothPadded[n:n+len(T)]
    Asmooth[0] = 0.
    return Asmooth

def ahe_from_tT(t,T,radius=50e-6,U=5e-5,Th=1e-5):
    '''
    t : numpy array
        time in sec
    T : numpy array
        temperature in Kelvin
    default radius and U from Guo et al 2021, GCA 310: 113-130.
'''
    Ma_to_s = 1e6 * 365.25 * 24 * 60 * 60
    t = t * Ma_to_s
    T = T + 273.15
    
    #print('radius, U, Th', radius, U, Th)
    
    ahe = helium_diffusion_models.calculate_he_age_meesters_dunai_2002(t, T, 
          radius, U, Th)
    
    # convert back to Ma and take only the final value
    return ahe[-1] / Ma_to_s


def AT_from_p_dA(p_dA,p_T,dT,Tsmooth):
    '''
    Generate arrays of age and temperature from input parameters.

    Parameters
    ----------
    p_dA : array-like
        Input parameters of age increments between each value of p_T.
    p_T : array-like
        Input parameters of temperature (deg C).
    dT : float
        Temperature increment in Temperature array.
    Tsmooth : float
        Temperature smoothing window width 
        (for sigma of Gaus smoothing window).

    Returns
    -------
    A : nd-array
        Age.
    T : nd-array
        Temperature.

    '''
    T = np.arange(p_T[0],p_T[-1]+dT,dT)
    p_A = p_A_from_p_dA(p_dA)
    A = A_from_p_A(p_A,T,dT,Tsmooth)
    return A,T

def aft_from_tT(t,T,obsArgs):
    '''
    Calculate age from time, temperature history

    Parameters
    ----------
    t : nd-array
        Time, forwards.
    T : nd-array
        Temperature at each time (usually getting less with time, cooling).
    obsArgs : tuple
        obs,obsErr,Dpar,p_T,dT,Tsmooth,Tdiff,geothermalGradient = obsArgs

    Returns
    -------
    result : nd-array
        Predicted AFT age, mean track length and standard deviation.
        [age, L_mean, L_sd].

    '''
    obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,geothermalGradient = obsArgs
    
    r = AFTannealingLib.simulate_AFT_annealing(t,T, 
#        kinetic_parameter='rmr0',kinetic_value=0.8, 
        kinetic_parameter='Dpar',kinetic_value=Dpar, 
#        kinetic_parameter='Dpar',kinetic_value=2.5, 
        method='Ketcham2007',
        initial_track_length = -9999, # Use appropriate default
#        initial_track_length = 16.3, # Durango is 16.3
#        use_fortran_algorithm=False,
        use_fortran_algorithm=True,
        )
    # return age, L_mean, L_sd
    return np.asarray(r[1:4])

def aft_from_p_dA(p_dA,obsArgs):    
    # returns individual values: aftAge, Lmean, Lstd
    
    t,T = tT_from_p_dA(p_dA,obsArgs)
    
    return aft_from_tT(t,T,obsArgs)

def aftGrad_from_p_dA(p_dA,obsArgs):
    '''
    Calculate AFT age, Lmean, Lstd, grad_aftAge, grad_Lmean, grad_Lstd
    Note that the grad_aftAge is the spatial age gradient, predicted assuming
    the geothermalGradient variable passed within obsArgs

    Parameters
    ----------
    p_dA : array-like
        Parameters of age differences.
    obsArgs : tuple
        obs,obsErr,Dpar,p_T,dT,Tsmooth,Tdiff,geothermalGradient = obsArgs

    Returns
    -------
    modelValues : nd-array
        Output model values.
        [AFT_age, Lmean, Lstd, grad_aftAge, grad_Lmean, grad_Lstd]
        Gradients are d/dT for changes in whole T(t) history

    '''
    obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,geothermalGradient = obsArgs
    
    t,T = tT_from_p_dA(p_dA,obsArgs)

    modelHighT = aft_from_tT(t,T + 0.5*Tdiff,obsArgs)
    modelLowT = aft_from_tT(t,T - 0.5*Tdiff,obsArgs)
    
    meanVals  = (modelHighT + modelLowT) / 2.0
    gradients = (modelHighT - modelLowT) / Tdiff
    
    modelValues = np.concatenate([meanVals,gradients])
    
    return modelValues

def convertTime(age,temperature):
    '''
    Time is forwards, but age is backwards, so flip and zero.
    Works for time -> age, or age -> time, and also does np.flip(T)
    
    Parameters
    ----------
    age : array-like
        Age or time.
    temperature : array-like
        Temperature.

    Returns
    -------
    time : array-like
        Age or time (opposite to input).
    temperature : array-like
        Temperature, flipped to match new time/age scale.
        
    '''
    return (age[-1] - np.flip(age)), np.flip(temperature)

def fitCoolingRate(p_dA,obsArgs):
    '''
    Generates age,T profile, finds model cooling rate at age of observations.
    Returns sum of weighted squared residuals

    Parameters
    ----------
    p_dA : nd-array
        Model adjustable parameters.
    coolRate : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,geothermalGradient = obsArgs
    dTdA_age = coolRate[:,0]
    dTdA_obs = coolRate[:,1]
    dTdA_err = coolRate[:,2]
    
    A,T = AT_from_p_dA(p_dA,p_T,dT,Tsmooth)
    Agrad = np.gradient(A)
    Tgrad = np.gradient(T)
    dTdA_mesh = Tgrad/Agrad
    dTdA_interpFunc = interp1d(A,dTdA_mesh,fill_value='extrapolate')
    dTdA_model = dTdA_interpFunc(dTdA_age)
    
    rss = np.sum(np.square((dTdA_obs-dTdA_model)/dTdA_err))
       
    return rss

def fitObs(p_dA,obsArgs,logfile):
    '''
    Fit observations to model defined by parameters [...dAi...].
    Return a single residual sum of squares scalar for optimisation

    Parameters
    ----------
    p_dA : nd-array
        dA parameters.
    obsArgs : tuple
        obs,obsErr,Dpar,p_T,dT,Tsmooth,Tdiff,geothermalGradient = obsArgs
        with:    
            obs : nd-array
                observations [age,Lmean,Lsd,ageGrad_m].
            obsErr : nd-array
                standard errors of observations.
            p_T : nd-array
                Temperatures parameters.
            dT : float
                Increment of temperature used in age model.
            Tsmooth : float
                Temperature window to smooth T(t) model.
            Tdiff : float
                Temperature window to calcualte gradient over.
            geothermalGradient : float
                dT/dz in degC/m, note z is up, so usually negative, e.g. -0.03

    Returns
    -------
    rss : float
        Residual sum of squared residuals sum(square((obs-model)/err)).

    '''
    obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,geothermalGradient = obsArgs
    
    modelVals = aftGrad_from_p_dA(p_dA,obsArgs)[0:4]
    # Convert age gradient prediction from Myr/C to Myr/m 
    # dA/dz = dA/dT * dT/dz
    modelVals[3] = modelVals[3] * geothermalGradient
    
    Q1 = np.sum(np.square((obs - modelVals)/obsErr))
    
    Q2 = fitCoolingRate(p_dA,obsArgs)
    
    rss = Q1 + Q2 * alpha/len(coolRate)
    
    #print(rss,modelVals)

    logfile.write(str(rss) + ',' + str(Q1) + ',' + str(Q2) + ',' + str(alpha) +
                  ',' + str(obs[0]) + ',' +
                  ','.join(map(str, p_dA)) + '\n')
    
    return rss

def fitParameters(p_dA,obsArgs,pBounds,filename='logfile.txt'):
    
    with open(filename,'a') as logfile:
    
        simplexMultipliers = simplexMultiplier_noUniformChange(p_dA) 
        
        p = np.array(p_dA)
        
        for m in simplexMultipliers:
    
            simplex = simplexNormalizeAge(m * p, obsArgs)
        
            result = minimize(fitObs,p_dA,args=(obsArgs,logfile),
                              bounds=pBounds,
                              method='Nelder-Mead',
                              options={'initial_simplex':simplex},
                              tol=0.2)
            
            p = result.final_simplex[0].mean(axis=0)
            qFit = result.fun
            
        return qFit,p
   
def load_model_output_bins():
    '''
    Loads the optimized model output for later analysis.
    Usage:
        obsArgs,model_output_bins = load_model_output_bins()

    Returns
    -------
    obsArgs : tuple 
        obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,geothermalGradient = obsArgs

    model_output_bins : tuple 
        aft_model_bins, p_dA_bins = model_output_bins
    '''
    # load best model parameters
    # Get summary parameters for help in making T(A) curves
    # model results
    log = pd.read_csv('logAll.csv')
    results = pickle.load(open('logOverviewResults.pickle','rb'))
    obsArgsList, pFit, obsPred = results
    obsArgs = obsArgsList[0]
    obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,geothermalGradient = obsArgs

    # Extract best models
    alphaChoice = 2
    p_dA_list = list()
    aft_model_list = list()
    for i in range(9):
        sel = (log['sample'] == i)
        logSel = log[sel & (log['alpha']==alphaChoice)]
        p_dA = logSel[logSel['rss']==logSel['rss'].min()].values[0,5:12]
        p_dA_list.append(p_dA)
        aft_model_list.append(aft_from_p_dA(p_dA,obsArgs)[0])
    aft_model_bins = np.asarray(aft_model_list)
    p_dA_bins = np.asarray(p_dA_list)
    model_output_bins = (aft_model_bins, p_dA_bins)
    return obsArgs,model_output_bins

def logCreate(p_dA,filename='logfile.txt'):
    with open(filename,'w') as logfile:
        logfile.write('rss,Q1,Q2,alpha,age')
        for i in range(len(p_dA)):
            logfile.write(',dA'+str(i))
        logfile.write('\n')
   
def p_A_from_p_dA(p_dA):
    '''
    Converts parameters from [...dAn] to [...Am], m = n+1
    '''
    return np.concatenate(([0.0],np.cumsum(p_dA)))

def p_dA_from_aft(age,obsArgs,model_output_bins):
    '''
    First run:
        obsArgs,model_output_bins = load_model_output_bins()
    '''
    obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,geothermalGradient = obsArgs
    aft_model_bins, p_dA_bins = model_output_bins
    aft_min = aft_model_bins[0]
    aft_max= aft_model_bins[-1]
    interp_p_dA_from_age = interp1d(aft_model_bins,p_dA_bins, axis=0)

    if age < aft_min:
        p_dA = interp_p_dA_from_age(aft_min)
        p_dA[1] -= (aft_min - age)
        
    elif age < aft_max : 
        p_dA = interp_p_dA_from_age(age)
    
    else:
        p_dA = interp_p_dA_from_age(aft_max)
        p_dA[2] += (age - aft_max)
        
    # check model AFT age is correct
    aft_model = aft_from_p_dA(p_dA,obsArgs)[0]
    dA = age - aft_model 
    if (p_dA[0] + dA) > 0 :
        p_dA[0] += dA
    else:
        p_dA[1] += dA
            
    return p_dA

def plot_p_dA(p_dA,ax,color='lightgrey'):
    A,T = AT_from_p_dA(p_dA,p_T,dT,Tsmooth)
    ax.plot(A,T,c=color)
    
def plotLogfile(filename='logfile.txt'):
    
    log = pd.read_csv('logfile.txt')
    #logSel = log[log['rss'] < 1.5*qFit] # may need F-test if can't fit
    log95 = log[log['Q1'] < 9.49] # chisq(4) @ 0.05 confidence
    log50 = log[log['Q1'] < 3.36] # chisq(4) @ 0.5 confidence
    p95 = log95.iloc[:,3:10].values
    p50 = log50.iloc[:,3:10].values
    pBest = log[log['rss']==log['rss'].min()].values[0,3:10]
    print(pBest)
    
    # T array is the same for all models, find A values for best parameters
    Abest,T = AT_from_p_dA(pBest,p_T,dT,Tsmooth)             
    
    fig, ax = plt.subplots(1)
    
    if len(p95) > 0:    
        A95 = list()
        for i in range(len(p95)):
            A,T = AT_from_p_dA(p95[i],p_T,dT,Tsmooth)
            A95.append(A)
        A95 = np.asarray(A95)   
        A95min = A95.min(axis=0)
        A95max = A95.max(axis=0)
        ax.fill_betweenx(T,A95min,A95max,color='pink')
    
    if len(p50) > 0:
        A50 = list()
        for i in range(len(p50)):
            A,T = AT_from_p_dA(p50[i],p_T,dT,Tsmooth)
            A50.append(A)
        A50 = np.asarray(A50)   
        A50min = A50.min(axis=0)
        A50max = A50.max(axis=0)
        ax.fill_betweenx(T,A50min,A50max,color='steelblue')
        
    #for A in Apaths: ax.plot(A,T,c='lightgrey')
    
    #ax.plot(Amin,T,c='pink')
    #ax.plot(Amax,T,c='pink')
    ax.plot(Abest,T,c='black')
            
    ax.set(xlabel='Age (Ma)',ylabel=r'Temperature ( $\degree$C)')
    ax.set(xlim=[0,150],ylim=[-10,160])
#    ax.set_xlim(0,150)
#    ax.set_ylim(-10,130)
    ax.invert_xaxis()
    ax.invert_yaxis()
    plt.tight_layout()
    plt.savefig('logfile.jpg')

def simplexMultiplier_noUniformChange(params,increase=1.5,decrease=0.5):
    '''
    Build an array of possible starting simplexes of parameter multipliers.
    
    This implementation excludes the possibility of a uniform increase to all 
    parameters, or a uniform decrease. See simplexMultiplier()
    
    Implementation example for pGuess array with len(pGuess)=nParams:
        simplexes = simplexMultiplier_noUniformChange(nParams) * pGuess

    Parameters
    ----------
    nParams : int
        Number of parameters.
    increase : float, optional
        Factor to increase parameter guess. The default is 1.25.
    decrease : float, optional
        Factor to decrease parameter guess. The default is 0.75.

    Returns
    -------
    simplexMultipliers : nd-array
        An array of simplexes. 
        Each simplex is an array of parameter multipliers.

    '''
    nParams = len(params)
    # Each simplex has one point with input guess.
    # We don't want all increase or all decrease so remove possibility: - 1
    # For other situations, we may want 
    # one = np.ones(nParams)
    one = np.ones(nParams-1)
    possibilities = list(one*increase) + list(one*decrease)
    # produce a list of possible p_dA multipliers for each new guess
    choices = list(set(itertools.permutations(possibilities,nParams)))

    simplexMultiplier = list()
    for i in range(0,len(choices),nParams):
        
        # build simplexes with one choice equal to last guess
        nodeList = [np.ones(nParams)]
        
        for j in range(nParams):
            # if we go past end of possible choices, start at beginning again
            nodeList.append(choices[np.mod(i+j, len(choices))])
                                    
        simplexMultiplier.append(nodeList)
    
    return np.asarray(simplexMultiplier)

def simplexMultiplier(params,increase=1.5,decrease=0.5):
    '''
    Build an array of possible starting simplexes of parameter multipliers.
    
    For implementation that excludes possibility of a uniform increase to all 
    parameters or a uniform decrease, see simplexMultiplier_noUniformChange()
    
    Implementation example for pGuess array with len(pGuess)=nParams:
        simplexes = pGuess * simplexMultiplier(nParams)

    Parameters
    ----------
    nParams : int
        Number of parameters.
    increase : float, optional
        Factor to increase parameter guess. The default is 1.25.
    decrease : float, optional
        Factor to decrease parameter guess. The default is 0.75.

    Returns
    -------
    simplexMultipliers : nd-array
        An array of simplexes. 
        Each simplex is an array of parameter multipliers.

    '''
    nParams = len(params)

    one = np.ones(nParams)
    possibilities = list(one*increase) + list(one*decrease)
    # produce a list of possible p_dA multipliers for each new guess
    choices = list(set(itertools.permutations(possibilities,nParams)))

    simplexMultiplier = list()
    for i in range(0,len(choices),nParams):
        
        # build simplexes with one choice equal to last guess
        nodeList = [np.ones(nParams)]
        
        for j in range(nParams):
            # if we go past end of possible choices, start at beginning again
            nodeList.append(choices[np.mod(i+j, len(choices))])
                                    
        simplexMultiplier.append(nodeList)
    
    return np.asarray(simplexMultiplier)
   
def simplexNormalizeAge(simplex,obsArgs): 
    
    obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,geothermalGradient = obsArgs
    ageTarget = obs[0]

    # first point already optimized so miss out
    for j in range(1,len(simplex)):
        ageModel = aft_from_p_dA(simplex[j],obsArgs)[0]
        simplex[j] = simplex[j] * ageTarget / ageModel           
    return simplex  

def tT_from_p_dA(p_dA,obsArgs):
    obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,geothermalGradient = obsArgs
    A,T = AT_from_p_dA(p_dA,p_T,dT,Tsmooth)
    return convertTime(A,T)
    
if __name__ == '__main__':
    
    #-------------------------------------------------------------------------
    tStart = time.process_time()
    
    #-------------------------------------------------------------------------
    # Array of observations [obs0,obs1...] where obs = [age,Lmean,Lsd,ageGrad_m]
    # Predicted ageGrad is ageGrad_T in Myr/C (negative sign)
    # To convert to Myr/m (spatial, positive)
    # dA/dz = dA/dT * dT/dz
    # But z positive up, T positive down; assume geothermalGradient = -0.03 C/m
    # ageGrad_m = - ageGrad_T * geothermalGradient
    # exhumation at 1 mm/yr = 1000 m/Myr gives dA/dz = 0.001 Myr/m
    obs =    np.array([50., 14.0, 1.7, 0.02 ]) #[age,Lmean,Lsd,ageGrad_m]
    obsErr = np.array([2.,   0.1, 0.1, 0.01])
    Dpar = 1.7
    coolRate = np.array([[40.,4.,2.],[70,1,0.3]])
    alpha = 1.
    obsArgs = (obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,geothermalGradient)
    # obs,obsErr,Dpar,p_T,dT,Tsmooth,Tdiff,geothermalGradient = obsArgs
    
    
    # Guess first parameters, uniform cooling
    # Ages at reference temperatures p_T
    p_A = obs[0] * p_T / (Tclosure - p_T[0])
    # Age differences - our choice of adjustable parameter
    p_dA = np.diff(p_A)
    # Bounding values of dA to control optimization
    pBounds = np.column_stack([0.01*p_dA,100.*p_dA])
    
    logCreate(p_dA)  
    
    qFit,pFit = fitParameters(p_dA,obsArgs,pBounds)
    
    #qFit,pFit = fitParameters(pFit,obsArgs)
    
    plotLogfile()
    
    print('time',time.process_time() - tStart)
    
    
