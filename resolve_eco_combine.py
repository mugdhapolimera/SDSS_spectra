# -*- coding: utf-8 -*-
"""
Created on Thu Mar 07 17:07:28 2019

@author: mugdhapolimera
"""

import pandas as pd
import os

#os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra')

resolve = pd.read_csv('RESOLVE_snr5.csv')
eco = pd.read_csv('ECO_snr5.csv')
eco = eco.rename(columns = {"NAME": "name"})

full = eco.copy()
notineco = (resolve['econame'] == 'notineco')
full = full.append(resolve[notineco])
full.index = full.name
print(full)        
#full.to_pickle('ECO+RESOLVE_filter_new.pkl')
full.to_csv('ECO+RESOLVE_snr5.csv')
print (len(full))