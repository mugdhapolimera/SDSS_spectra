# -*- coding: utf-8 -*-
"""
Created on Wed Mar 06 13:30:23 2019

@author: mugdhapolimera


wget -nv -r -nH --cut-dirs=7 -i speclist.txt -B https://data.sdss.org/sas/dr12/sdss/spectro/redux/26/spectra/lite/
"""

import pandas as pd

sdss = pd.read_csv('SDSS_spectral_plateid.csv')

speclist  = open('speclist.txt', 'w')

for plate, mjd, fiberid in zip(sdss['plate'],sdss['mjd'],sdss['fiberid']):
    count+=1    
    speclist.write("{}/spec-{}-{}-{}.fits \n".format(plate, plate, mjd, fiberid))

speclist.close()