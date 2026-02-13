#!/usr/bin/env python
# python3
"""
    Analyse aft tracks
"""
__author__ = "Rupert Sutherland"

import numpy as np
import pandas as pd
import pyproj as proj
from scipy.spatial import KDTree
from scipy.stats import binned_statistic

def grad(position,a):
    # linear regression of scalar values y = y0 + Gx = Mx
    # so design matrix M has rows 1 + x + y + z
    M = np.column_stack([np.ones_like(position[:,0]),position])
    degreesOfFreedom = M.shape[0] - M.shape[1]
    MtM_inverse = np.linalg.inv((M.T).dot(M))
    solutionMatrix = MtM_inverse.dot(M.T)
    pHat = solutionMatrix.dot(a)
    aHat = M.dot(pHat)
    aResidual = (a - aHat)[:,np.newaxis] # add dimension for matrix operations
    aS2 = (aResidual.T).dot(aResidual)/degreesOfFreedom
    pCov = aS2 * MtM_inverse
    gradVec = pHat[1:]
    gradVecCov = pCov[1:,1:]
    gradMag = np.sqrt(gradVec.dot(gradVec))
    gradDirection = gradVec/gradMag
    gradMagErr = np.sqrt((gradDirection.T).dot(gradVecCov).dot(gradDirection))
    return gradMag,gradMagErr,gradVec


wdir = '../data/'
file1 = wdir + 'aft_tam.csv'

aft = pd.read_csv(file1)

sample = aft['sample'].values
lon = aft['lon'].values
lat = aft['lat'].values
a = aft['age'].values
s = aft['age_serr'].values
z = aft['z'].values
nData = len(a)

# Project lon/lat into x,y in metres
# EPSG 4326: lon,lat WGS84
# EPSG 32758: UTM 58S central longitude 165E
crs = proj.CRS.from_epsg(32758)
transformer = proj.Transformer.from_crs("EPSG:4326","EPSG:32758",always_xy=True)
x,y = transformer.transform(aft['lon'],aft['lat'])
position = np.column_stack([x,y,z])

# The things to find
idList = list()
lonList = list()
latList = list()
nList = list()
xyzList = list()
aMeanList = list()
aMeanErrList = list()
aGradDirXList = list()
aGradDirYList = list()
aGradDirZList = list()
aGradMagList = list()
aGradMagErrList = list()
# also then compute and append to output: exRate exRateErr

# use KDTree to step through each point and find nearest neighbours
# then compute means and gradients from those points
# nLocal is the nukmber of points used to compute means and gradients
nMin = 6 
nMax = 17
nRange = np.arange(nMin,nMax+1,1)

for nLocal in nRange:
    
    tree = KDTree(position)
    
    for i in range(nData):
        
        # k is array of indices of nearest points
        distances,k = tree.query(position[i],k=nLocal)
        
        gradMag,gradMagErr,gradVec = grad(position[k],a[k])
        
        gradDir = gradVec/gradMag
        
        idList.append(sample[i])
        lonList.append(lon[i])
        latList.append(lat[i])
        nList.append(nLocal)
        #xyzList.append(position[k].mean(axis=0).tolist())
        xyzList.append(position[k].mean(axis=0))
        aMeanList.append(a[k].mean())
        aMeanErrList.append(a[k].std(ddof=1)/np.sqrt(nLocal))
        aGradMagList.append(gradMag)
        aGradMagErrList.append(gradMagErr)
        aGradDirXList.append(gradDir[0])
        aGradDirYList.append(gradDir[1])
        aGradDirZList.append(gradDir[2])

df = pd.DataFrame({'sample':idList,
                   'lon':lonList,
                   'lat':latList,
                   'n':nList,
                   'xyzMean':xyzList,
                   'ageMean':aMeanList,
                   'ageMeanErr':aMeanErrList,
                   'ageGrad':aGradMagList,
                   'ageGradErr':aGradMagErrList,
                   'aGradDirX':aGradDirXList,
                   'aGradDirY':aGradDirYList,
                   'aGradDirZ':aGradDirZList   })   

df['exRate'] = np.reciprocal(df['ageGrad'])

# u = 1/x , du/dx = -1/x**2 , du = - dx / x**2
df['exRateErr'] = df['ageGradErr'] / df['ageGrad']**2

df.to_csv(wdir + 'aftGrad.csv', index=False)

#---------------------------------------------------
# Select data then bin by age and create summary output file

# Restrict to NVL
sel = df['lat'] > -76

print('ageMean min,max ',df['ageMean'][sel].min(),',',df['ageMean'][sel].max())

# Quality control results
sel = sel & (df['exRateErr'] > 0.)
sel = sel & (df['exRateErr'] < 1000.)
sel = sel & (df['exRateErr'] < (1. * df['exRate']))

# Don't trust a result unless age gradient is sub-vertical upwards
sel = sel & (df['aGradDirZ'] > 0.9)  

selGood = sel 

# Now look at exhumation rate binned by time (Myr)
binSize = 5.
binAgeMax = 100.

aBins = np.arange(int(df['ageMean'][sel].min()), binAgeMax, binSize)   
binAgeMedian = aBins[:-1] + 0.5 * np.diff(aBins)

# build new dataframe for age-binned data

# Each value of nLocal is needed for each age bin
binN = np.repeat(nRange,len(binAgeMedian))
# We want a value for each age bin for each value of nLocal in nRange
binAge = np.tile(binAgeMedian,len(nRange))

# Add columns for mean, standard error, nValuesUsed
nBinned = list()
ageMean = list()
ageMeanErr = list()
exRate = list()
exRateErr = list()

for nLocal in nRange:

    sel = selGood & (df['n'] == nLocal)
    
    n = binned_statistic(df['ageMean'][sel], df['ageMean'][sel],
                         statistic='count',bins=aBins)[0]
    
    nBinned.append(n)  
    
    ageMean.append(binned_statistic(df['ageMean'][sel],
                                    df['ageMean'][sel],
                                    statistic='mean',
                                    bins=aBins)[0]
                   )
    ageMeanErr.append(binned_statistic(df['ageMean'][sel],
                                    df['ageMean'][sel],
                                    statistic='std',
                                    bins=aBins)[0] /
                  np.sqrt(n-1)                  
                   )
    exRate.append(binned_statistic(df['ageMean'][sel],
                                    df['exRate'][sel],
                                    statistic='mean',
                                    bins=aBins)[0]
                   )
    exRateErr.append(binned_statistic(df['ageMean'][sel],
                                    df['exRate'][sel],
                                    statistic='std',
                                    bins=aBins)[0] /
                  np.sqrt(n-1)                  
                   )

dfBin = pd.DataFrame({'nLocal' : binN, 
                      'ageBin' : binAge,
                      'nBinned': np.concatenate(nBinned),
                      'ageMean': np.concatenate(ageMean),
                      'ageMeanErr': np.concatenate(ageMeanErr),
                      'exRate': np.concatenate(exRate),
                      'exRateErr': np.concatenate(exRateErr)
                      })

dfBin['ageGrad'] = np.reciprocal(dfBin['exRate'])

# u = 1/x , du/dx = -1/x**2 , du = - dx / x**2
dfBin['ageGradErr']  =  dfBin['exRateErr'] * dfBin['ageGrad']**2

dfBin.to_csv(wdir + 'aftGrad_nvl_binMean.csv',index=False)

#---------------------------------------------------
# Create a summary of average values from all nLocal
# add a summary and give it an id of n=0

ageMean = list()
ageMeanErr = list()
ageGrad = list()
ageGradErr = list()
exRate = list()
exRateErr = list()

for a in binAgeMedian:
    sel = (dfBin['ageBin'] == a)
    ageMean.append(np.nanmean(dfBin[sel]['ageMean'].values))
    # sample std
    ageMeanErr.append(np.nanstd(dfBin[sel]['ageMean'].values,ddof=1))
    # mean of age variances 
    #ageMeanErr.append(np.sqrt(np.mean(np.square(dfBin[sel]['ageMeanErr'].values))))
    exRate.append(np.nanmean(dfBin[sel]['exRate'].values))
    # find sample standard deviation
    exRateErr.append(np.nanstd(dfBin[sel]['exRate'].values,ddof=1))

dfAv = pd.DataFrame({'nMin': np.ones_like(binAgeMedian)*nMin,
                     'nMax': np.ones_like(binAgeMedian)*nMax,
                     'ageBin': binAgeMedian,
                     'ageMean': ageMean,
                     'ageMeanErr': ageMeanErr,
                     'exRate': exRate,
                     'exRateErr': exRateErr
                     })

dfAv['ageGrad'] = np.reciprocal(dfAv['exRate'])
dfAv['ageGradErr']  =  dfAv['exRateErr'] * dfAv['ageGrad']**2

dfAv.to_csv(wdir + 'aftGrad_nvl_binSummary.csv',index=False)
