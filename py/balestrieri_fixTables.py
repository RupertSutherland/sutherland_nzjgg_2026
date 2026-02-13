#!/usr/bin/env python
# python3
"""
    Balestrieri et al 2020
    Geoscience Frontiers 11 (2020) 1841–1858
    Tables reformat
"""
__author__ = "Rupert Sutherland"

import numpy as np
from scipy.stats import chi2
import pandas as pd


wdir = '../data/'

infile1 = wdir + 'aft_balestrieriTable1.txt'
outfile1 = wdir + 'aft_balestrieriTable1_fixed.csv'
infile2 = wdir + 'ahe_balestrieriTable2.txt'
outfile2 = wdir + 'ahe_balestrieriTable2_fixed.csv'
outfile3 = wdir + 'ahe_balestrieriTable2_meanage.csv'

source = 'Balestrieri et al. (2020)'

def decDeg(x):
    y = x.split('_')
    degs = float(y[0])
    mins = float(y[1][0:2])
    secs = float(y[1][3:])
    return round(degs + mins/60. + secs/3600., 5)

def list2str(myList):
    return str(myList).replace(',','').replace('[','').replace(']','')

# Fix Table 1 AFT data
with open(infile1,'r') as f:
    lines = f.readlines()
with open(outfile1,'w') as fout:
    fout.write(
               'lon' + ',' + 
               'lat' + ',' + 
               'z' + ',' + 
               'age' + ',' + 
               'age_serr' + ',' + 
               'Lm' + ',' + 
               'Lm_serr' + ',' + 
               'L_sd' + ',' + 
              # 'prob' + ',' +
              # 'n_grains' + ',' +
               'Dpar' + ',' +
               'source' + ',' +
               'sample' + '\n')
    
    for line in lines:
        word = line.strip().split(' ')
        sample = word[0].replace('_','-')
        lat = str(-decDeg(word[3]))
        lon = str(decDeg(word[4]))
        z = word[5]
        n_grains = word[12]
        prob = word[13].replace('<','')
        age = word[14]
        age_serr = word[16]
        Lm = word[18]
        Lm_serr = word[20]
        L_sd = word[21]
        if Lm == '0' :
            Lm = '0'
            Lm_serr = '999'
            L_sd = '999'
            Dpar = word[22].split('(')[0]
        else:
            Dpar = word[23].split('(')[0]
        
        fout.write( 
                   lon + ',' + 
                   lat + ',' + 
                   z + ',' + 
                   age + ',' + 
                   age_serr + ',' + 
                   Lm + ',' + 
                   Lm_serr + ',' + 
                   L_sd + ',' + 
               #    prob + ',' +
               #    n_grains + ',' + 
                   Dpar + ',' + 
                   source + ',' + 
                   sample +'\n')
    
df1 = pd.read_csv(outfile1)

# Fix Table 2 AHE data
with open(infile2,'r') as f:
    lines = f.readlines()
    # add EOF string
    lines = lines + ['EOF']
    
with open(outfile2,'w') as f2out:
    with open(outfile3,'w') as f3out:
        f2out.write(      'lon' + 
                    ',' + 'lat' + 
                    ',' + 'z' + 
                    ',' + 'age' + 
                    ',' + 'age_serr' + 
                    ',' + 'include' + 
                    ',' + 'source' + 
                    ',' + 'sample' + 
                    '\n')   
        f3out.write(      'lon' + 
                    ',' + 'lat' + 
                    ',' + 'z' + 
                    ',' + 'age' + 
                    ',' + 'age_serr' + 
                    ',' + 'age_ssd' + 
                    ',' + 'n' + 
                    ',' + 'source' + 
                    ',' + 'sample' + 
                    ',' + 'chiSq'  + 
                    ',' + 'p_chiSq' + 
                    ',' + 'ages' + 
                    ',' + 'ages_errs' + 
                    ',' + 'includes' + '\n')

        sample = 'no samples yet'
        ageList = list()
        errList = list()
        includeList = list()
        
        for line in lines:
            word = line.strip().split(' ')
            sample1 = word[0]
            sample2 = word[0][:-1]
            if sample1 != 'EOF':
                lon = str(df1[df1['sample']==sample2]['lon'].values[0])
                lat = str(df1[df1['sample']==sample2]['lat'].values[0])
                z =   str(df1[df1['sample']==sample2]['z'].values[0])
                age = word[11]
                age_serr = word[13].replace('a','')
                include = str(int(not('a' in word[13])))
            f2out.write(      lon + 
                        ',' + lat + 
                        ',' + z + 
                        ',' + age + 
                        ',' + age_serr + 
                        ',' + include + 
                        ',' + source + 
                        ',' + sample1 +
                        '\n')   
            
            # calculate pooled sample age after reading all grain ages
            if sample2 == sample:
                ageList.append(float(age))
                errList.append(float(age_serr))
                includeList.append(bool(int(include)))
                
            else:
                if sample != 'no samples yet':
                    #print(includeList)
                    i = np.asarray(includeList,dtype=bool)
                    ages = np.asarray(ageList,dtype=float)[i]
                    ageMean = np.mean(ages)
                    n = len(ages)
                    ages_ssd = np.std(ages,ddof=1)
                    ageMean_serr = ages_ssd/np.sqrt(n)
                    ages_errs = np.asarray(errList,dtype=float)[i]
                    w = ages_errs**-2
                    ageWeightedMean = np.sum(np.multiply(w,ages))/np.sum(w)
                    chisq = np.sum(np.multiply(w,np.square(ages-ageMean)))
                    p_chisq = chi2(n).cdf(chisq)
                    f3out.write(      lon + 
                                ',' + lat + 
                                ',' + z + 
                                ',' + str(np.round(ageMean,1)) + 
                                ',' + str(np.round(ageMean_serr,1)) + 
                                ',' + str(np.round(ages_ssd,1)) + 
                                ',' + str(n) + 
                                ',' + source + 
                                ',' + sample + 
                                ',' + str(np.round(chisq,1))  + 
                                ',' + str(p_chisq) + 
                                ',' + list2str(ageList) + 
                                ',' + list2str(errList) + 
                                ',' + list2str(includeList) + '\n') 
                                  
                sample = sample2
                ageList = [float(age)]
                errList = [float(age_serr)]            
                includeList = [bool(include)]
        



