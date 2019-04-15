# -*- coding: utf-8 -*-
"""
Created on Thu Mar 07 17:07:28 2019

@author: mugdhapolimera
"""

import pandas as pd
import os

os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra')

resolve = pd.read_pickle('RESOLVE_filter_new.pkl')
eco = pd.read_pickle('ECO_filter_new.pkl')
eco = eco.rename(columns = {"name": "NAME"})

full = eco.copy()
notineco = (resolve['econame'] == 'notineco')
full = full.append(resolve[notineco])
        
full.to_pickle('ECO+RESOLVE_filter_new.pkl')
full.to_csv('ECO+RESOLVE_filter_new.csv')