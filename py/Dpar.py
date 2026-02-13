#!/usr/bin/env python
# python3
"""
    Plot AFT Dpar data.

"""
__author__ = "Rupert Sutherland"

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import rc
rc("pdf", fonttype=42)
rc("font",size=10,family='Arial')
rc("figure",figsize=(4, 4))

wdir = '../data/'

aft = pd.read_csv(wdir + 'aft_tam.csv')

aft = aft[aft['lat'] > -76]
aft = aft[aft['Dpar'] > 0.]

plt.hist(aft['Dpar'],cumulative=-1,color='lightgrey',edgecolor='black',
         bins=np.arange(1.2,3.,0.2))

meanDpar = aft['Dpar'].mean()

plt.axvline(x=1.7,label='mean Dpar 1.7',color='black')

plt.xticks(np.arange(1.2,3.,0.4))

plt.xlabel('Dpar ($\mu$m)')
plt.ylabel('Number of measurements')

plt.legend()

plt.tight_layout()

plt.savefig(wdir+'Dpar.pdf')

print(aft['Dpar'].mean())
