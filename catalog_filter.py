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

eco = 0
resolve = 1
portsmouth = 0
jhu = 1
if eco: 
    inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_full_blend_dext_new.pkl'
    outputfile = "C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_full_snr5_port"
    df = pd.read_pickle(inputfile)
    print (len(df))
    mgas = df.logmgas
    mstars = df.logmstar
    mbary = 10**mgas + 10**mstars
    inobssample = (((df.grpcz >= 3000.) & (df.grpcz <= 7000.))) #& \
    #((np.log10(mbary) > 9.2))) #(df.absrmag < -17.3) & 


#In RESOLVE Sample filtering
if resolve:
    inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_blend_dext_new.pkl'
    outputfile = "C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_snr5_port"
    df = pd.read_pickle(inputfile)
    print len(df)
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
    inobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & \
    (((flinsample | (np.log10(mbary) > 9.0)) & infall) | \
            ((flinsample | (np.log10(mbary) > 9.2)) & inspring))
#    inobssample = (((grpcz >= 4500.) & (grpcz <= 7000.)) & \
#    (np.log10(mbary) > 9.2))


#for ECO
#Line FLuxes Filtering
floor = 10**-3
ceil = 1e5
df = df[inobssample]
print len(df)
if portsmouth:
    nii = df['Flux_NII_6583']
    nii_sum = (df['Flux_NII_6583']+ df['Flux_NII_6547'])*3./4
    nii_sum_err = (np.sqrt(df['Flux_NII_6547_Err']**2 + df['Flux_NII_6583_Err']**2))*3./4
    oiii = df['Flux_OIII_5006']
    oiii_err = df['Flux_OIII_5006_Err']
    h_alpha = df['Flux_Ha_6562']
    h_alpha_err = df['Flux_Ha_6562_Err']
    h_beta = df['Flux_Hb_4861']
    h_beta_err = df['Flux_Hb_4861_Err']
    oi = df['Flux_OI_6300']
    oi_err = df['Flux_OI_6300_Err']
    sii_sum = df['Flux_SII_6716'] + df['Flux_SII_6730']
    sii_sum_err = np.sqrt(df['Flux_SII_6716_Err']**2 + df['Flux_SII_6730_Err']**2)

if jhu:
    nii = df['nii_6584_flux']
    if 'nii_6548_flux' in df.keys():
        nii_sum = (df['nii_6584_flux']+ df['nii_6548_flux'])*3./4
        nii_sum_err = (np.sqrt(df['nii_6584_flux_err']**2 + df['nii_6548_flux_err']**2))*3./4
    else:
        nii_sum = df['nii_6584_flux']
        nii_sum_err = df['nii_6584_flux_err']
    # note that the ratio uses only the stronger line, but for S/N reasons we add
    # the weaker and multiply by 3/4 since Chris Richardson says the canonical
    # line ratio is 3:1 (this needs to be updated with a more precise number)
    oiii = df['oiii_5007_flux']
    oiii_err = df['oiii_5007_flux_err']
    h_alpha = df['h_alpha_flux']
    h_alpha_err = df['h_alpha_flux_err']
    h_beta = df['h_beta_flux']
    h_beta_err = df['h_beta_flux_err']
    oi = df['oi_6300_flux']
    oi_err = df['oi_6300_flux_err']
    if 'sii_6717_flux' in df.keys():
        sii_sum = df['sii_6717_flux'] + df['sii_6731_flux']
    
        sii_sum_err = np.sqrt(df['sii_6717_flux_err']**2 + df['sii_6731_flux_err']**2)
    else:
        sii_sum = df['sii_6731_flux']
    
        sii_sum_err = df['sii_6731_flux_err']

gooddata = ((h_alpha > floor) & (nii_sum > floor) & (oiii > floor) & (oi > floor) &
            (sii_sum > floor) & (h_beta > floor)  & (h_alpha_err > floor) & 
            (nii_sum_err > floor) & (oiii_err > floor) & 
            (oi_err > floor) & (sii_sum_err > floor) & 
            (h_alpha < ceil) & (nii_sum < ceil) & (oiii < ceil) & (oi < ceil) &
            (sii_sum < ceil) & (h_beta < ceil) & 
            ~np.isnan(h_alpha) & ~np.isnan(nii_sum) & ~np.isnan(oiii) & 
            ~np.isnan(oi) & ~np.isnan(sii_sum) & ~np.isnan(h_beta))

df = df[gooddata]
cut = 5
snr = ((h_alpha > cut*h_alpha_err) & (nii_sum > cut*nii_sum_err) & 
       (oiii > cut*oiii_err) & (oi > cut*oi_err) & (sii_sum > cut*sii_sum_err) 
       & (h_beta > cut*h_beta_err))

print len(df)
df = df[snr]
print len(df) 
#flags = pd.read_csv('C:/Users/mugdhapolimera/github/BPT/resolve_emlineclass_bpt1_new.csv')
#flags.index = flags['galname']

#df = df[flags['defstarform']]
#df.to_pickle(outputfile+".pkl")
#df.to_csv(outputfile+".csv")
#print len(df)

#old
#df = df[~np.isnan(df.h_alpha_flux) & (df.h_alpha_flux > floor)
#        & (df.h_alpha_flux < ceil)]
#df = df[~np.isnan(df.oiii_5007_flux) & (df.oiii_5007_flux > floor) 
#        & (df.oiii_5007_flux < ceil)]
#df = df[~np.isnan(df.nii_6584_flux) & (df.nii_6584_flux > floor)
#        & (df.nii_6584_flux < ceil)]
#df = df[~np.isnan(df.nii_6548_flux) & (df.nii_6548_flux > floor)
#        & (df.nii_6548_flux < ceil)]
#df = df[~np.isnan(df.h_beta_flux) & (df.h_beta_flux > floor)
#        & (df.h_beta_flux < ceil)]
#df = df[~np.isnan(df.oi_6300_flux) & (df.oi_6300_flux > floor)
#       & (df.oi_6300_flux < ceil)]
#df = df[~np.isnan(df.sii_6717_flux) & (df.sii_6717_flux > floor)
#        & (df.sii_6717_flux < ceil)]
#df = df[~np.isnan(df.sii_6731_flux) & (df.sii_6731_flux > floor)
#        & (df.sii_6731_flux < ceil)]
##df = df[~np.isnan(df.heii_4685_flux_port_ext) & 
##        (df.heii_4685_flux_port_ext > floor) & (df.heii_4685_flux_port_ext < ceil)]
#snr = 5
#df = df[(df.h_beta_flux/df.h_beta_flux_err >= snr)]
#df = df[df.h_alpha_flux/df.h_alpha_flux_err >= snr]
#df = df[df.oiii_5007_flux/df.oiii_5007_flux_err >=snr ]
#df = df[df.nii_6584_flux/df.nii_6584_flux_err >= snr ]
#df = df[df.nii_6548_flux/df.nii_6548_flux_err >= snr ]
#df = df[df.sii_6717_flux/df.sii_6717_flux_err >=snr ]
#df = df[df.sii_6731_flux/df.sii_6731_flux_err >=snr ]
##print len(df)
#df = df[df.oi_6300_flux/df.oi_6300_flux_err >= snr ]
