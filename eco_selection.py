# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 12:18:43 2019

@author: mugdhapolimera
"""

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
import matplotlib.pyplot as plt
eco = 1

inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_blend_dext_new.pkl'
df = pd.read_pickle(inputfile)
print (len(df))
mgas = df.logmgas
mstars = df.logmstar
mbary = 10**mgas + 10**mstars
#inobssample = ((df.grpcz >= 3000.) & (df.grpcz <= 7000.)) & ((df.absrmag < -17.3) & (np.log10(mbary) > 9.2))
inobssample = ((df.grpcz >= 4500.) & (df.grpcz <= 7000.) & (np.log10(mbary) > 9.2))
#inobssample = ((df.grpcz >= 3000.) & (df.grpcz <= 7000.)) & ((np.log10(mbary) > 9.2))

plt.figure()
plt.plot(df.absrmag, np.log10(mbary),'b.')
total = np.sum(np.log10(mbary) > 9.0)
df = df[inobssample]
mgas = df.logmgas
mstars = df.logmstar
mbary = 10**mgas + 10**mstars
sel = float(np.sum(np.log10(mbary) > 9.0))
print sel/total
plt.plot(df.absrmag, np.log10(mbary),'r.')
plt.plot(-17.3+0*np.linspace(8,12), np.linspace(8,12),'k')
plt.plot(np.linspace(-24,-16),9.2+0*np.linspace(-24,-16),'k')
plt.xlim(-16,-24)
plt.ylim(8,12)

#for ECO

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
df = df[df.nii_6548_flux/df.nii_6548_flux_err >= 3 ]
df = df[df.oi_6300_flux/df.oi_6300_flux_err >= 2 ]
df = df[df.sii_6717_flux/df.sii_6717_flux_err >=2 ]
df = df[df.sii_6731_flux/df.sii_6731_flux_err >=2 ]
print len(df)
#flags = pd.read_csv('C:/Users/mugdhapolimera/github/BPT/resolve_emlineclass_bpt1_new.csv')
#flags.index = flags['galname']

#df = df[flags['defstarform']]
#df.to_pickle("C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_filter_new.pkl")
#df.to_csv("C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_filter_new.csv")
#print len(df)
