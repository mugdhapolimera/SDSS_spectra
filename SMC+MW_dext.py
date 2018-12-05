# -*- coding: utf-8 -*-
"""
Created on Sun Dec 02 16:30:50 2018

SMC + Milky Way Blended Extinction Correction

This code blends SMC and MW based extinction correction, based on the 
galaxy mass. 

    - SMC extinction curve (Gordon03) for dwarf galaxies with M < 10^9
    - MW extinction curve (O’Donnell94) for giant galaxies with M > 10^10
    - A linear combination of SMC and MW extinction for galaxies with 
      10^9 < M < 10^10
                A = logM % 9 , B = 1-A
                Dereddened flux = A * SMC extinction + B * MW extinction

@author: mugdhapolimera
"""
from astropy.table import Table
import numpy as np
import pandas as pd
smcfile = 'C:\Users\mugdhapolimera\github\SDSS_Spectra\RESOLVE_SDSS_full_smcdext.fits'
smc0 = Table.read(smcfile, format='fits')
smc = smc0.to_pandas()
smc.index = smc.NAME
mwfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_SDSS_full_dext.fits'
mw0 = Table.read(mwfile, format='fits')
mw = mw0.to_pandas()
mw.index = mw.NAME
resolve = pd.read_csv('RESOLVE_liveOctober2018.csv')
resolve.index = resolve.name
ext_corr = smc.copy()

for gal in smc.index.values:
    if resolve.logmstar.loc[gal] < 9:
        ext_corr.loc[gal] = smc.loc[gal]
    elif resolve.logmstar.loc[gal] > 10:
        ext_corr.loc[gal] = mw.loc[gal]
    else:
        A = resolve.logmstar.loc[gal] % 9
        B = 1 - A
        print 'Blending for Galaxy ', gal, A, B
        print smc.loc[gal][1:73]
        ext_corr.loc[gal][1:73] = A * smc.loc[gal][1:73] + B * mw.loc[gal][1:73]
print ext_corr