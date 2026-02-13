#!/usr/bin/env python
# python3
"""
    Amalgamate new aft/ahe data
"""
__author__ = "Rupert Sutherland"

import numpy as np
import pandas as pd


wdir = '../data/'

infile1 = wdir + 'allApatites.xlsx'
infile2 = wdir + 'aft_balestrieriTable1_fixed.csv'
infile3 = wdir + 'ahe_balestrieriTable2_fixed.csv'

outfile1 = wdir + 'aft_tam.csv'
outfile2 = wdir + 'ahe_tam.csv'
'''
lonMin = 166.0
lonMax = 172.0
latMin = -72.8
latMax = -70.0
'''
lonMin = 150.0
lonMax = 180.0
latMin = -90.0
latMax = -60.0

df1 = pd.read_excel(infile1,sheet_name='allapatites',header=1)
df1['L_sd'] = df1['track_length_2sig']/2.0
# make an 'include' column with integer=1 values
df1['include'] = (df1['Top_Depth/Sample_Elevation'] -
                  df1['Top_Depth/Sample_Elevation'] + 1 )

dfs1 = df1[(df1['Latitude']>latMin)  & (df1['Latitude']<latMax) &
           (df1['Longitude']>lonMin) & (df1['Longitude']<lonMax)]

#-------------------------------------------------

# exclude any data not clearly AFT ages with track lengths
dfs1aft = dfs1[dfs1['track_length_mean']>0]

# Identify columns that will match columns in infile2 (Balestrieri Table1)
dfs1aftColumns = ['Longitude', 'Latitude','Top_Depth/Sample_Elevation',
              'Age_Ma','1sigma','track_length_mean','track_length_std_err',
              'L_sd','Dpar','Laboratory','Collection_Concat']

# Now load infile2 (Balestrieri Table1) and concatenate with other data
df2aft = pd.read_csv(infile2)
aftCols = df2aft.columns

df_aft = pd.DataFrame(dfs1aft[dfs1aftColumns].values,columns=aftCols)

df_aft = pd.concat([df_aft,df2aft])

df_aft.to_csv(outfile1,index=False)

#-------------------------------------------------

dfs1ahe = dfs1[dfs1['Geochron_Method']=='U-Th-He']

dfs1aheColumns = ['Longitude', 'Latitude','Top_Depth/Sample_Elevation',
              'Age_Ma','1sigma','include','Laboratory','Collection_Concat']

df3ahe = pd.read_csv(infile3)
aheCols = df3ahe.columns

df_ahe = pd.DataFrame(dfs1ahe[dfs1aheColumns].values,columns=aheCols)

df_ahe = pd.concat([df_ahe,df3ahe])

df_ahe.to_csv(outfile2,index=False)
#-------------------------------------------------


