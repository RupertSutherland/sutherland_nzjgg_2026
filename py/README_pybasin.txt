Search for modifications by Rupert to catch errors in pybasin/lib/AFTannealingLib.py 
Comments start with '##RS'

Change made to function simulate_AFT_annealing()

        ##RS - catch Nan values and change to 0.0
        rmp[np.isnan(rmp)] = 0.0

