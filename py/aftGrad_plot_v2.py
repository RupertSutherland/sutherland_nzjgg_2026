#!/usr/bin/env python
# python3
"""
    Analyse aft tracks
"""
__author__ = "Rupert Sutherland"

import numpy as np
import pandas as pd
from scipy.stats import binned_statistic

import matplotlib.pyplot as plt
from matplotlib import rc
rc("pdf", fonttype=42)
rc("font",size=7,family='Arial')
rc("figure",figsize=(4,3),dpi=300)

wdir = '../data/'

df = pd.read_csv(wdir + 'aftGrad_nvl_binMean.csv')
df2 = pd.read_csv(wdir + 'aftGrad_nvl_binSummary.csv')

nMin = 6 
nMax = 17
nRange = np.arange(nMin,nMax+1,1)

# start plotting
xmin = 25 
xmax = 90
plt.xlim(xmin,xmax)
plt.ylim(0,800)
plt.xticks(range(xmin,xmax,5)) 
# c20y
plt.plot([43.5,43.5],[400,800],color='black',marker=None,
         linewidth=1,linestyle='--') 
# c24y
plt.plot([52.9,52.9],[400,800],color='black',marker=None,
         linewidth=1,linestyle='--') 
plt.xlabel('Age (Ma)')
plt.ylabel('Exhumation rate (m/Ma)')


cmap = plt.cm.get_cmap('coolwarm')
col = cmap(np.linspace(0,1,len(nRange)))

for nLocal in nRange:
    
    sel = (df['nLocal'] == nLocal) & (df['exRate'] > 0.)

    '''
    plt.figure(1)
    plt.errorbar(aMean[sel],exRate[sel],xerr=aStd[sel],yerr=exRateErr[sel],
                 fmt='.',ecolor='lightgrey',color='black',markersize=3)    
    plt.figure(2)
    plt.errorbar(aBinned,exRateBinned,yerr=sBinned,
                 fmt='.',ecolor='pink',color='darkred')    
    '''
    #plt.plot(aBinned,exRateBinned,color='steelblue',marker=None,linewidth=1)    
    #plt.plot(aBinned,exRateBinned + sBinned,color='lightgrey',marker=None)    
    #plt.plot(aBinned,exRateBinned - sBinned,color='lightgrey',marker=None) 
    
    '''
    plt.scatter(df['ageMean'][sel], df['exRate'][sel],
             marker='o', s=9,
             c=col[nLocal-nMin]) 
    plt.errorbar(df['ageMean'][sel], df['exRate'][sel],
                 xerr=df['ageMeanErr'][sel],yerr=df['exRateErr'][sel],
                 fmt='.',ecolor='lightgrey',color='grey',markersize=1)    
    '''
    plt.errorbar(df['ageMean'][sel], df['exRate'][sel],
                 xerr=df['ageMeanErr'][sel],yerr=df['exRateErr'][sel],
                 fmt='.',markersize=4,
                 ecolor=col[nLocal-nMin],color=col[nLocal-nMin])    
 

plt.plot(df2['ageMean'], df2['exRate'],
         color='black',marker=None,linewidth=1)    

plt.plot(df2['ageMean'], df2['exRate'] + df2['exRateErr'],
         color='black',linestyle='--',marker=None,linewidth=1) 
   
plt.plot(df2['ageMean'], df2['exRate'] - df2['exRateErr'],
         color='black',linestyle='--',marker=None,linewidth=1) 
   
'''
plt.errorbar(df2['ageMean'], df2['exRate'],
             xerr=df2['ageMeanErr'],yerr=df2['exRateErr'],
             fmt='.',ecolor='steelblue',color='steelblue',markersize=1)    
'''
#plt.plot(aMean,exMean,color='black',marker=None,linewidth=2)
plt.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(nMin,nMax)),
             label='N points', ticks=np.arange(8,17,2),
             orientation='vertical', shrink=0.3)

plt.savefig(wdir + 'xRateNVL.jpg',bbox_inches='tight')
plt.savefig(wdir + 'xRateNVL.pdf',bbox_inches='tight')
