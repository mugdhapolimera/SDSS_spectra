
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 07 17:07:28 2019

@author: mugdhapolimera
"""

import pandas as pd
import os
import numpy as np

os.chdir('C:/Users/mugdhapolimera/github/SDSS_spectra/mid_ir/')

resolve = pd.read_csv('RESOLVE_WISE_good.csv')
#mgas = resolve.logmgas
#mstars = resolve.logmstar
#mbary = 10**mgas + 10**mstars
#print(np.sum(mbary > 10**9.2))
eco = pd.read_csv('ECO_WISE_good.csv')#_xray_chandra_new ECO/SEL/ECO_full_snr5_dext_nsa.csv')
#mgas = eco.logmgas
#mstars = eco.logmstar
#mbary = 10**mgas + 10**mstars
#print(np.sum(mbary > 10**9.2))
eco = eco.rename(columns = {"NAME": "name_x"})
 
full = resolve.copy()
notinresolve = (eco['resname'] == 'notinresolve')
full = full.append(eco[notinresolve])
full.index = full.name
#mgas = full.logmgas
#mstars = full.logmstar
#mbary = 10**mgas + 10**mstars
#full = full[mbary > 10**9.2]
#print(full)        
#full.to_pickle('ECO+RESOLVE_filter_new.pkl')
full.to_csv('ECO+RESOLVE_WISE_good.csv')# _xray_chandra_new ECO+RESOLVE_snr5_dext_nsa.csv')
print (len(full))