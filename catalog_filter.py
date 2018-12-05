# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 11:46:29 2018

@author: mugdhapolimera

Filter all sources in the RESOLVE catalog according to the following criteria:

In Observation Sample:
    - Group recession velocity between 4500 and 7000 km/s
    - In Sample Flag = 1
    - Log Baryonic Mass
        - Fall: Greater than 9
        - Spring: Greater than 9.2
    - JHU Flag- reliable = 1 

Line fluxes and Errors:
    - All lines used for bpt (SELs)- not NaN
    - Greater than 0, less than 10000 (arbitrary max flux value)
    - H_beta : at least 3 sigma
    - All other lines : at least 2 sigma

"""
import numpy as np
import pandas as pd 

inputfile = 'C:/Users/mugdhapolimera/github/izi/RESOLVE_SDSS_full.pkl'
df = pd.read_pickle(inputfile)
print (len(df))
#In RESOLVE Sample filtering
ra=df.radeg
dec=df.dedeg
flinsample = df.fl_insample
grpcz = df.grpcz
cz = df.cz
infall = (ra > 22*15.) | (ra < 3*15.)
inspring = (ra > 8.75*15.) & (ra < 15.75*15.)
mgas = df.logmgas
mstars = df.logmstar
mbary = 10**mgas + 10**mstars
inobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & (((flinsample | (np.log10(mbary) > 9.0)) & infall) | ((flinsample | (np.log10(mbary) > 9.2)) & inspring))


#Line FLuxes Filtering
floor = 0#10**-2
ceil = 1e5
df = df[inobssample]

df = df[~np.isnan(df.h_alpha_flux_ext) & (df.h_alpha_flux_ext > floor)
        & (df.h_alpha_flux_ext < ceil)]
df = df[~np.isnan(df.oiii_5007_flux_ext) & (df.oiii_5007_flux_ext > floor) 
        & (df.oiii_5007_flux_ext < ceil)]
df = df[~np.isnan(df.nii_6584_flux_ext) & (df.nii_6584_flux_ext > floor)
        & (df.nii_6584_flux_ext < ceil)]
df = df[~np.isnan(df.nii_6548_flux_ext) & (df.nii_6548_flux_ext > floor)
        & (df.nii_6548_flux_ext < ceil)]
df = df[~np.isnan(df.h_beta_flux_ext) & (df.h_beta_flux_ext > floor)
        & (df.h_beta_flux_ext < ceil)]
#df = df[~np.isnan(df.oi_6300_flux_ext) & (df.oi_6300_flux_ext > floor)
#        & (df.oi_6300_flux_ext < ceil)]
#df = df[~np.isnan(df.sii_6717_flux_ext) & (df.sii_6717_flux_ext > floor)
#        & (df.sii_6717_flux_ext < ceil)]
#df = df[~np.isnan(df.sii_6731_flux_ext) & (df.sii_6731_flux_ext > floor)
#        & (df.sii_6731_flux_ext < ceil)]
print len(df)
df.to_pickle("C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_SDSS_inobssample.pkl")

#df = df[~np.isnan(df.heii_4685_flux_port_ext) & 
#        (df.heii_4685_flux_port_ext > floor) & (df.heii_4685_flux_port_ext < ceil)]

df = df[(df.h_beta_flux_ext/df.h_beta_flux_ext_err >= 3)]
df = df[df.h_alpha_flux_ext/df.h_alpha_flux_ext_err >= 2]
df = df[df.oiii_5007_flux_ext/df.oiii_5007_flux_ext_err >=2 ]
df = df[df.nii_6584_flux_ext/df.nii_6584_flux_ext_err >= 2 ]
df = df[df.nii_6548_flux_ext/df.nii_6548_flux_ext_err >= 2 ]
df = df[df.oi_6300_flux_ext/df.oi_6300_flux_ext_err >= 2 ]
df = df[df.sii_6717_flux_ext/df.sii_6717_flux_ext_err >=2 ]
df = df[df.sii_6731_flux_ext/df.sii_6731_flux_ext_err >=2 ]

#df.to_pickle("C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_SDSS_filtered.pkl")

