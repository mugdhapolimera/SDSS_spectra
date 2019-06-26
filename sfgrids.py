# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 14:37:46 2019

@author: mugdhapolimera
"""

import pandas as pd
import os

os.chdir("C:/Users/mugdhapolimera/github/izi/")
#gridfile = r'Richardson-0-0_1-0agn-BPASS-Binary-CSF-n=1e2-40.0Myr-NichollsCE.csv'

#gridfile = 'Richardson-0-0_1-0agn-BPASS-Binary-CSF-n=1e2-40.0Myr-NichollsCE-D_G-RR14_Fstar_0_3.csv'
#gridfile = 'Richardson-0-0_1-0agn-BPASS-Binary-SSP-n=1e2-1.0Myr-NichollsCE-D_G-RR14_Fstar_0_3.csv'
gridfile = 'Richardson-0-0_1-0agn-BPASS-Binary-SSP-n=1e2-20Myr-NichollsCE-D_G-RR14_Fstar_0_3.csv'

grid = pd.read_csv(gridfile)
grid = grid[grid["AGNFRAC"] == 0]
grid = grid[(grid["LOGQ"] > 6.9) & (grid["LOGQ"] < 8.9)]

grid.to_csv('Richardson-0-sf-BPASS-Binary-SSP-n=1e2-20Myr-NichollsCE-D_G-RR14_Fstar_0_3.csv')

