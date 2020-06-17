# -*- coding: utf-8 -*-
"""
Created on Thu Mar 07 17:07:28 2019

@author: mugdhapolimera
"""

import pandas as pd
import os

os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra')

resolve = pd.read_csv('RESOLVE_snr5_master_bary.csv')
eco = pd.read_csv('ECO_snr5_master_bary.csv')
eco = eco.rename(columns = {"NAME": "name"})

full = resolve.copy()
notinresolve = (eco['resname'] == 'notinresolve')
full = full.append(eco[notinresolve])
full.index = full.name
print(full)        
#full.to_pickle('ECO+RESOLVE_filter_new.pkl')
full.to_csv('ECO+RESOLVE_snr5_master_bary.csv')
print (len(full))