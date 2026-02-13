#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Process logfiles to output results in useful form
pickle dump of results = (obsArgsList,pFit,obsPred)
To restore later:
    import pickle
    results = pickle.load(open('logOverviewResults.pickle','rb'))
    obsArgsList, pFit, obsPred = results

and synthesis of all possible paths, confidence regions

@author: rupert
"""

import numpy as np
import pandas as pd

import pickle

import aftModel_v2 as model

def strToArray(s):
    '''
    Converts a string output onto one line with str(array) i.e. of form
    '[[1 2][3 4]]' and returns np.array([1,2],[3,4])
    

    Parameters
    ----------
    s : str
        format '[[1 2][3 4]]'.
        Can be '[[1 2][3 4]]\n'

    Returns
    -------
    a : nd-array
        format np.array([1,2],[3,4]).

    '''
    from ast import literal_eval
    sSingleSpaced = ' '.join(s.split())
    sFormatted = sSingleSpaced.replace('[ ' , '[')
    sFormatted = sFormatted.replace(' ]' , ']')
    sFormatted = sFormatted.replace(' ' , ',')
    sFormatted = sFormatted.replace('][' , '],[')
    #print(sFormatted)
    return np.array(literal_eval(sFormatted))

#a = strToArray('[[1 2][3 4]]\n')

wdir = '../data/'


'''
Each block of logOverview.txt has format...
NEW_MODEL
obs
obsErr
Dpar
coolRate
alpha
p_T
dT
Tsmooth
Tdiff
geothermalGradient
qFit 0.0002880614479274299
pFit [22.94428836  5.2639607   0.08368587  0.08341437  0.0663763  10.80819383 21.06328781]
time(s): 117.42861342299966

obsArgs =  (obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,geothermalGradient)
obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,geothermalGradient = obsArgs
'''

obsArgsList = list()
p = list()
obsPred = list()

with open('logOverview.txt','r') as f:
    lines = f.readlines()

for i in range(len(lines)):
    if lines[i] == 'NEW_MODEL\n':
        obs = strToArray(lines[i+1])
        obsErr = strToArray(lines[i+2])
        Dpar = float(lines[i+3])
        coolRate = strToArray(lines[i+4])
        alpha = float(lines[i+5])
        p_T = strToArray(lines[i+6]) 
        dT = float(lines[i+7])
        Tsmooth = float(lines[i+8])
        Tdiff = float(lines[i+9])
        geothermalGradient = float(lines[i+10])
        pFit = strToArray(lines[i+12].replace('pFit ',''))
        
        obsArgs = (obs,obsErr,Dpar,coolRate,alpha,p_T,dT,Tsmooth,Tdiff,
                        geothermalGradient)
        
        pred = model.aftGrad_from_p_dA(pFit,obsArgs)[0:4]
        
        obsArgsList.append(obsArgs)
        p.append(pFit)
        obsPred.append(pred)
        
p = np.asarray(p)
obsPred = np.asarray(obsPred)

results = (obsArgsList,p,obsPred)

pickle.dump(results,open('logOverviewResults.pickle','wb'))

#--------------------------------------------------------------------------
# Now get all individual log files and process for each datum, to find range
# of possible parameters
   
for alpha in range(10):
    
    for sample in range(11):
        
        filename = 'logfile' + str(alpha) + '{:02d}'.format(sample) + '.txt'
            
        df = pd.read_csv(filename)       
        #df['alpha'] = alpha
        df['sample'] = sample
        
        if alpha==0 and sample==0:
            dfAll = df
        else:
            dfAll = pd.concat([dfAll,df])
    
dfAll.to_csv('logAll.csv',index=False)    
