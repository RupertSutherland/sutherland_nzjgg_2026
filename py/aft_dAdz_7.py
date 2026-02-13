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

import matplotlib.pyplot as plt
from matplotlib import rc
rc("pdf", fonttype=42)
rc("font",size=7,family='Arial')
rc("figure",figsize=(4,3),dpi=300)

def grad(position,a):
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
aft = aft[aft['lat'] > -76]

sample = aft['sample'].values
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
nList = list()
xyzList = list()
aMeanList = list()
aMeanErrList = list()
aGradDirList = list()
aGradMagList = list()
aGradMagErrList = list()
# also then compute and append to output: exRate exRateErr

# use KDTree to step through each point and find nearest neighbours
# then compute means and gradients from those points
nMin = 5 
nMax = 20
nRange = np.arange(nMin,nMax+1,1)

for nLocal in nRange:
    
    tree = KDTree(position)
    
    for i in range(nData):
        
        # k is array of indices of nearest points
        distances,k = tree.query(position[i],k=nLocal)
        
        gradMag,gradMagErr,gradVec = grad(position[k],a[k])
        
        idList.append(sample[i])
        nList.append(nLocal)
        xyzList.append(position[k].mean(axis=0).tolist())
        aMeanList.append(a[k].mean())
        aMeanErrList.append(a[k].std(ddof=1)/np.sqrt(nLocal))
        aGradMagList.append(gradMag)
        aGradMagErrList.append(gradMagErr)
        aGradDirList.append((gradVec/gradMag).tolist())

names = ['sample','n','xyz','ageMean','ageMeanErr','ageGrad','ageGradErr',
         'ageGradVec']

df = pd.DataFrame({'sample':idList,
                   'n':nList,
                   'xyz':xyzList,
                   'ageMean':aMeanList,
                   'ageMeanErr':aMeanErrList,
                   'ageGrad':aGradMagList,
                   'ageGradErr':aGradMagErrList,
                   'ageGradDir':aGradDirList })   

df['exRate'] = np.reciprocal(df['ageGrad'])

# u = 1/x , du/dx = -1/x**2 , du = - dx / x**2
df['exRateErr'] = df['ageGradErr'] / df['ageGrad']**2

df.to_csv(wdir + 'aftGrad.csv', index=False)

