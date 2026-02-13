Codes to go with publication 
Eocene tectonic change in the South Pacific caused exhumation of the Transantarctic Mountains in northern Victoria Land, Antarctica
R. Sutherland, P.G. Fitzgerald, N. Mortimer.
New Zealand Journal of Geology and Geophysics 2026

---------------------------------
aft_dAdz.py DEPRICATED: aftGrad.py

Calculate spatial age gradient (hence exhumation rate).

Input: aft_tam.csv
Output: aftGrad.csv

Uses scipy.spatial.KDTree to find neigbours.
Uses linear regression to find gradient.

---------------------------------
aft_dAdz_map.py DEPRICATED: aftGrad.py

Plots graph of exhumation rate as a function of time.

Also does calculations for specific n.
Input: aft_tam.csv
Output: xRate*.jpg



---------------------------------
aft_tracks.py

Bins data to calculate summary AFT stats, track lengths, etc.
 
Input: aft_tam.csv
Output: aft_tracks_nvl.jpg

---------------------------------
aft-ahe_data.py

Creates data file of AFT and AHE ages on same samples.
 
Input: aft_tam.csv, ahe_tam.csv
Output: aft-ahe_data.csv

---------------------------------
aft-ahe.py

Analyse differences in AFT and AHE ages.
 
Input: aft-ahe_data.csv

---------------------------------
AFTannealingLib.py

python library for apatite fission track annealing algorithms
from pyBasin - Elco Luijendijk, Goettingen University

---------------------------------
aftDataMeans.py

Load and plot summary AFT data. 
Output file has summary stats for each age bin (primary observations to model).
 
Input: aftGrad_nvl_binSummary.csv, aft_nvl_tracks.csv
Output: aft_allData.csv

---------------------------------
aftGrad.py

Calculates spatial age gradient of AFT data, binned by age.
 
Input: aft_tam.csv
Output: aftGrad.csv, aftGrad_nvl_binMean.csv, aftGrad_nvl_binSummary.csv

---------------------------------
aftGrad_plot.py

Calculates spatial age gradient of AFT data, binned by age.
 
Input: aftGrad_nvl_binMean.csv, aftGrad_nvl_binSummary.csv
Output: xRateNVL.jpg, xRateNVL.pdf

---------------------------------
aftModel.py

Functions to model cooling histories to fit AFT data.
Uses AFTannealingLib from pyBasin.
 
---------------------------------
aftModelAll.py

Model cooling histories to fit AFT data for each age bin [i].
Uses aftModel.py
 
Input: aft_allData.csv
Output: logfile[alpha][age_bin].txt, logOverview.txt

---------------------------------
amalgamateData.py

Amalgamate aft/ahe data.
 
Input: allApatites.xlsx, aft_balestrieriTable1_fixed.csv, ahe_balestrieriTable2_fixed.csv
Output: aft_tam.csv, ahe_tam.csv

---------------------------------
balestrieri_fixTables.py

Balestrieri et al 2020; Geoscience Frontiers 11 (2020) 1841–1858
Tables reformat 

Input: aft_balestrieriTable1.txt, ahe_balestrieriTable2.txt
Output: aft_balestrieriTable1_fixed.csv, ahe_balestrieriTable2_fixed.csv, ahe_balestrieriTable2_meanage.csv

---------------------------------
Dpar.py

Plot AFT Dpar data.
Fig S1

Input: aft_tam.csv
Output: Dpar.pdf

---------------------------------
logfileVariance.py

variance-covariance estimated from rss estimations

Input: logfile.txt
Output: 

---------------------------------
plot_allData.py
NOTE - Old version, see below.

load and plot summary AFT data to model (3 subplots).

Input: aft_allData.csv
Output: aftSummary.jpg

---------------------------------
plotAFT_allData.py

load and plot summary NVL AFT data to model (3 subplots).
Fig 4.

Input: aft_tam.csv
Output: aft_allData.jpg

---------------------------------
plotOverview.py
NOTE - Old version, see below.

load and plot overview log.

Input: logOverview.txt
Output: aftOverview_allData.jpg

---------------------------------
plotOverview_allData.py

load and plot overview log.

Input: aft_allData.csv, logOverviewResults.pickle 
Output: aftOverview_allData.jpg

---------------------------------
plotOverview_AT.py

Plot overview log: age-temperature. 

Input: logOverviewResults.pickle
Output: aftOverview_AT.jpg

---------------------------------
plotOverview_confidence.py

Plot overview log: age-temperature; with confidence region and different alpha. 

Input: logOverviewResults.pickle
Output: aftOverview_samples.jpg

---------------------------------
plotOverviewPreprocess.py

Process logfiles to output results in useful form
pickle dump of results = (obsArgsList,pFit,obsPred)

Input: logOverview.txt 
Output: logOverviewResults.pickle

---------------------------------

