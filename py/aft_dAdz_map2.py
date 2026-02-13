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
a = aft['age'].values
s = aft['age_serr'].values
z = aft['z'].values
nData = len(a)
# Project lon/lat into x,y in metres
# EPSG 4326: lon,lat WGS84
# EPSG 32758: UTM 58S central longitude 165E
#crs = proj.CRS.from_epsg(32758)
transform_lonlat_utm58S = proj.Transformer.from_crs(
    "EPSG:4326","EPSG:32758",always_xy=True)
transform_utm58S_lonlat = proj.Transformer.from_crs(
    "EPSG:32758","EPSG:4326",always_xy=True)
x,y = transform_lonlat_utm58S.transform(aft['lon'],aft['lat'])

position = np.column_stack([x,y,z])

# use KDTree to step through each point and find nearest neighbours
nLocal = 7
tree = KDTree(position)
aMean = np.array([])
aStd = np.array([])
g = np.array([])
gErr = np.array([])
xMean = list()
for i in range(nData):
    distances,k = tree.query(position[i],k=nLocal)
    gradMag,gradMagErr,gradVec = grad(position[k],a[k])
    aMean = np.append(aMean,[a[k].mean()])
    aStd = np.append(aStd,[a[k].std(ddof=1)/np.sqrt(nLocal)])    
    g = np.append(g,[gradMag])
    gErr = np.append(gErr,[gradMagErr])
    xMean.append(position[k].mean(axis=0))

xMean = np.asarray(xMean)
exRate = np.reciprocal(g)
# u = 1/x , du/dx = -1/x**2 , du = - dx / x**2
exRateErr = gErr / g**2

lon,lat = transform_utm58S_lonlat.transform(xMean[:,0],xMean[:,1])

df = pd.DataFrame(data=xMean,columns=['x','y','z',])
df['lon'] = lon
df['lat'] = lat
df['exhumeRate'] = exRate
df['exRateErr'] = exRateErr
df['ageMean'] = aMean
df['ageMeanSerr'] = aStd
df.to_csv(wdir + 'exhumeRateNVL.csv')

sel = (exRateErr < 10000)
#---------------
aBins = np.arange(int(aMean.min()),aMean.max()+5,3)   
aBinned = binned_statistic(aMean[sel],aMean[sel],
                           statistic='mean',bins=aBins)[0]
exRateBinned = binned_statistic(aMean[sel],exRate[sel],
                                statistic='mean',bins=aBins)[0]
varMeanBinned =  binned_statistic(aMean[sel],exRateErr[sel]**2,
                                  statistic='mean',bins=aBins)[0]
nMeanBinned =  binned_statistic(aMean[sel],exRateErr[sel],
                                statistic='count',bins=aBins)[0]    
sBinned = np.sqrt(varMeanBinned/nMeanBinned)
'''
plt.figure(1)
plt.errorbar(aMean[sel],exRate[sel],xerr=aStd[sel],yerr=exRateErr[sel],
             fmt='.',ecolor='lightgrey',color='black',markersize=3)    
plt.figure(2)
plt.errorbar(aBinned,exRateBinned,yerr=sBinned,
             fmt='.',ecolor='pink',color='darkred')    
'''
plt.plot(aBinned,exRateBinned,color='darkred',marker=None,linewidth=1)    
#plt.plot(aBinned,exRateBinned + sBinned,color='lightgrey',marker=None)    
#plt.plot(aBinned,exRateBinned - sBinned,color='lightgrey',marker=None)    

#plt.savefig(wdir + 'xRateMapNVL.jpg',bbox_inches='tight')
