#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot all aft data

@author: rupert
"""

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import rc
rc("pdf", fonttype=42)
rc("font",size=6,family='Arial')
rc("figure",figsize=(3, 4))

wdir = '../data/'
file1 = wdir + 'aft_tam.csv'

aft = pd.read_csv(file1)

aft = aft[aft['Lm'] > 0.]
aft = aft[aft['lat'] > -76]

#------------------------------------------
df = pd.read_csv(wdir + 'aft_allData.csv')

fig, ax = plt.subplots(2)

ax[0].set(ylabel=r'Mean track length ( $\mu$m)')
ax[0].errorbar(aft['age'],aft['Lm'],
               xerr=aft['age_serr'],yerr=aft['Lm_serr'], 
               fmt='.',markersize=2,color='black',ls='none',
               ecolor='grey',elinewidth=1)

ax[1].set(ylabel=r'Track length s.d. ( $\mu$m)')
ax[1].errorbar(aft['age'],aft['L_sd'],
               xerr=aft['age_serr'], 
               fmt='.',markersize=2,color='black',ls='none',
               ecolor='grey',elinewidth=1)
ax[1].set(xlabel=r'Age (Ma)')

ax[0].set_xlim(0,120)
ax[1].set_xlim(0,120)
#ax[0].text(5,15,'A')
#ax[1].text(5,3.5,'B')
ax[0].set_title('A',loc='left')
ax[1].set_title('B',loc='left')


plt.tight_layout()

fig.savefig(wdir+'aft_allData.jpg',dpi=300)
fig.savefig(wdir+'aft_allData.pdf')

