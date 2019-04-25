# -*- coding: utf-8 -*-
"""
Created on Wed Mar 06 13:30:23 2019

@author: mugdhapolimera

Upload list of objects, ra, dec to https://skyserver.sdss.org/dr12/en/tools/crossid/crossid.aspx (Choose SpecObjAll table)

Use the following SQL snippet:

    SELECT s.specobjid, s.ra, s.dec, s.plate, s.mjd, s.fiberid
    FROM #upload u
          JOIN #x x ON x.up_id = u.up_id
          JOIN SpecObjAll s ON s.specObjID = x.specObjID

Download the resulting table as SDSS_spectral_plateid.csv

Create speclist.txt as below and run the following command from terminal
wget -nv -r -nH --cut-dirs=7 -i speclist.txt -B https://data.sdss.org/sas/dr12/sdss/spectro/redux/26/spectra/lite/
"""

import pandas as pd

sdss = pd.read_csv('SDSS_spectral_plateid.csv')

speclist  = open('speclist.txt', 'w')

for plate, mjd, fiberid in zip(sdss['plate'],sdss['mjd'],sdss['fiberid']):
       
    speclist.write("%04d/spec-%04d-%d-%04d.fits \n" %(plate, plate, mjd, fiberid))

speclist.close()
