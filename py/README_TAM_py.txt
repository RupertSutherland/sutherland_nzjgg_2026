Northern Victoria Land (NVL) Transantarctic Mountains (TAM) analysis.

Code in TAM/py
Raw and processed data files in TAM/data
Raw input is allApatites.xlsx and aft_balestrieriTable1.txt and ahe_balestrieriTable2.txt

First organize data

balestrieri_fixTables.py
allApatites.xlsx and aft_balestrieriTable1.txt and ahe_balestrieriTable2.txt  
-> 
aft_balestrieriTable1_fixed.csv & ahe_balestrieriTable2_fixed.csv & ahe_balestrieriTable2_meanage.csv

amalgamateData4.py
infile1 = wdir + 'allApatites.xlsx'
infile2 = wdir + 'aft_balestrieriTable1_fixed.csv'
infile3 = wdir + 'ahe_balestrieriTable2_fixed.csv'
-> 
outfile1 = wdir + 'aft_tam.csv'
outfile2 = wdir + 'ahe_tam.csv'

aft_tracks1.py
Bin aft data by age to give summary track lengths & sd
file1 = wdir + 'aft_tam.csv'
file2 = wdir + 'ahe_tam.csv' 
-> 
aft_nvl_tracks.csv

aftGrad.py
Calculate spatial age gradients and reciprocals = exhumation rates
Project locations, then linear regression on data grouped using KDTree. nMin = 6 ; nMax = 17
aft_tam.csv
->
aftGrad.csv = all data, value at every datum, every possible number of nLocal values used in regression
aftGrad_nvl_binMean.csv = NVL data binned by age, every possible number of nLocal values used in regression
aftGrad_nvl_binSummary.csv = NVL data binned by age, sample std from each nLocal (6-17 values). Note uncertainty introduced by data geometry more than original value of uncertainty.

aftGrad_plot_v2.py 
Generate plot of age gradients for paper
aftGrad_nvl_binMean.csv
aftGrad_nvl_binSummary.csv
->
xRateNVL.pdf

aftDataMeans.py 
Resample all nvl track length data from 10 Myr means to 5 Myr, like ageGrad data, to create master dataset to model.
aft_nvl_tracks.csv
aftGrad_nvl_binSummary.csv
->
aft_allData.csv

plot_allData.py
Load and plot nvl data to model
aft_allData.csv
->
aftSummary.pdf

AFTannealingLib.py 
Kinetics from pybasin. Use function AFTannealingLib.simulate_AFT_annealing() 
Can use f2py to compile calculate_reduced_AFT_lengths.f90 to make it run faster. See separate notes.

aftMod02_fitData_v1.py 
Working test of fitting an individual observation. Best run on linux vm from within spyder, using compiled f90.
NOTE: aft_from_tT(t,T,obsArgs) calculates an age by calling AFTannealingLib.simulate_AFT_annealing(). Check options here!





