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

inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_full_blend_dext.pkl'
df = pd.read_pickle(inputfile)
print (len(df))
#In RESOLVE Sample filtering
#ra=df.radeg
#dec=df.dedeg
#flinsample = df.fl_insample
#grpcz = df.grpcz
#cz = df.cz
#infall = (ra > 22*15.) | (ra < 3*15.)
#inspring = (ra > 8.75*15.) & (ra < 15.75*15.)
#mgas = df.logmgas
#mstars = df.logmstar
#mbary = 10**mgas + 10**mstars
#inobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & (((flinsample | (np.log10(mbary) > 9.0)) & infall) | ((flinsample | (np.log10(mbary) > 9.2)) & inspring))
#for ECO
inobssample = (df.absrmagcorr < -17.3)

#Line FLuxes Filtering
floor = 10**-3
ceil = 1e5
df = df[inobssample]
print len(df)
df = df[~np.isnan(df.h_alpha_flux) & (df.h_alpha_flux > floor)
        & (df.h_alpha_flux < ceil)]
df = df[~np.isnan(df.oiii_5007_flux) & (df.oiii_5007_flux > floor) 
        & (df.oiii_5007_flux < ceil)]
df = df[~np.isnan(df.nii_6584_flux) & (df.nii_6584_flux > floor)
        & (df.nii_6584_flux < ceil)]
df = df[~np.isnan(df.nii_6548_flux) & (df.nii_6548_flux > floor)
        & (df.nii_6548_flux < ceil)]
df = df[~np.isnan(df.h_beta_flux) & (df.h_beta_flux > floor)
        & (df.h_beta_flux < ceil)]
df = df[~np.isnan(df.oi_6300_flux) & (df.oi_6300_flux > floor)
        & (df.oi_6300_flux < ceil)]
df = df[~np.isnan(df.sii_6717_flux) & (df.sii_6717_flux > floor)
        & (df.sii_6717_flux < ceil)]
df = df[~np.isnan(df.sii_6731_flux) & (df.sii_6731_flux > floor)
        & (df.sii_6731_flux < ceil)]
#df = df[~np.isnan(df.heii_4685_flux_port_ext) & 
#        (df.heii_4685_flux_port_ext > floor) & (df.heii_4685_flux_port_ext < ceil)]
print len(df)
df = df[(df.h_beta_flux/df.h_beta_flux_err >= 3)]
df = df[df.h_alpha_flux/df.h_alpha_flux_err >= 3]
df = df[df.oiii_5007_flux/df.oiii_5007_flux_err >=3 ]
df = df[df.nii_6584_flux/df.nii_6584_flux_err >= 3 ]
#df = df[df.nii_6548_flux/df.nii_6548_flux_err >= 2 ]
#df = df[df.oi_6300_flux/df.oi_6300_flux_err >= 2 ]
#df = df[df.sii_6717_flux/df.sii_6717_flux_err >=2 ]
#df = df[df.sii_6731_flux/df.sii_6731_flux_err >=2 ]
#print len(df)
#flags = pd.read_csv('C:/Users/mugdhapolimera/github/BPT/eco_emlineclass_bpt1.csv')
#flags.index = flags['galname']

#df = df[flags['defstarform']]
df.to_pickle("C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_filter.pkl")
print len(df)
