#!/usr/bin/env python
# python3
__author__ = "Rupert Sutherland"
"""
    variance-covariance estimated from rss estimations
"""
import numpy as np
import pandas as pd


dat = pd.read_csv('logfile.txt')

# Residual sum of weighted squared residuals (chi-square 4 d.f.)
q = dat['rss'].values
df = 4.
# Parameter values at corresponding q values
p = dat.iloc[:,1:].values

# best estimate from mean value of good fits (within 10% of qBest)
pBest = p[q < 1.1*q.min()].mean(axis=0)

# parameter residuals 
dp = p - pBest

# Find a Hessian matrix estimate from each residual  
# scale each c matrix according to the q value
# the q value must be divided by df, so normalized to an expected value of 1
H = list()
for i in range(len(dp)):
    # turn each dp[i] into a column vector
    v = np.column_stack(dp[i]).T
    # dpi*dpj components
    dpipjMat = np.dot(v,v.T)
    # hessian is second derivative d2q/(dpi*dpj)
    H.append(2.*q[i] * np.reciprocal(dpipjMat))
H = np.asarray(H)

# select values suitable for estimation of Hessian
# chisq(4) @ 0.01 significance, 99% confidence is 13.28
# chisq(4) @ 0.05 significance, 95% confidence is 9.49
# chisq(4) @ 0.90 significance is 1.06
# Use 0.01 to 0.90 signficance, 10-99% confidence
s = ((q > 1.06) & (q < 13.28))

hessian = np.mean(H[s] , axis=0)

# the hessian matrix is the inverse of the covariance matrix
covP = np.linalg.inv(hessian)

