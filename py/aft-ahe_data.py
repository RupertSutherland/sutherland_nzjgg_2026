#!/usr/bin/env python
# python3
"""
    Creates data file of AFT and AHE ages on same samples
"""
__author__ = "Rupert Sutherland"

import numpy as np
import pandas as pd

from scipy.stats import norm
from scipy.stats import poisson

import matplotlib.pyplot as plt
from matplotlib import rc
rc("pdf", fonttype=42)
rc("font",size=7,family='Arial')
rc("figure",figsize=(4,3))



wdir = '../data/'

file1 = wdir + 'aft_tam.csv'
file2 = wdir + 'ahe_tam.csv'

aft = pd.read_csv(file1)
ahe = pd.read_csv(file2)

# for each aft age, look for all ahe ages within 11 m (0.0001 deg lat)
# Same samples should have dx=dy=0 
maxDist = 0.0001 

with open(wdir+'aft-ahe_data.csv', 'w') as f:
    
    f.write('age_AFT,age_AHE,id_AFT,id_AHE,check\n')
    
    for i in range(len(aft)):
        for j in range(len(ahe)):
            # skip if more than maxDist apart
            if (abs(aft['lon'][i] - ahe['lon'][j]) < maxDist and
                abs(aft['lat'][i] - ahe['lat'][j]) < maxDist) : 
                
                # write to output file if very close together
                f.write('{},{},{},{},{}\n'.format(aft['age'][i], ahe['age'][j],
                    aft['sample'][i], ahe['sample'][j],
                    aft['sample'][i] in ahe['sample'][j]))
            
            


