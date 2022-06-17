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

def catfilter(eco, resolve, sdsscat, saveoutput, hasnr, bpt1snr, he2,
              inputfile, outputfile, cut, he2cut):
    #Filter the sample as per survey definition

    if eco: 
        df = pd.read_csv(inputfile)
        df.index = df.name
        print (len(df))
        mgas = df.logmgas
        mstars = df.logmstar
        mbary = 10**mgas + 10**mstars
        ineco = (130.05 < df.radeg) & (df.radeg < 237.45)
        inobssample = (((df.grpcz >= 3000.) & (df.grpcz <= 7000.)) & 
                       (np.log10(mbary) > 9.2) & ineco)#\
    
    if resolve:    
        df = pd.read_csv(inputfile)
        df.index = df.name
        print(len(df))
        ra=df.radeg
        dec=df.dedeg
        flinsample = df.fl_insample
        grpcz = df.grpcz
        cz = df.cz
        infall = ((ra > 22*15.) | (ra < 3*15.)) & ((dec > -1.25) & (dec < 1.25))
        inspring = (ra > 8.75*15.) & (ra < 15.75*15.) & ((dec > 0) & (dec < 5))
        mgas = df.logmgas
        mstars = df.logmstar
        mbary = 10**mgas + 10**mstars
        inobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & \
        (((flinsample | (np.log10(mbary) > 9.0)) & infall) | \
                ((flinsample | (np.log10(mbary) > 9.2)) & inspring))
#        inobssample = (((grpcz >= 4500.) & (grpcz <= 7000.)) & \
#        (np.log10(mbary) > 9.2))
        df[inobssample].to_csv("RESOLVE_inobssample.csv")
    
    
    #Line FLuxes Filtering
    floor = 10**-3
    ceil = 1e5
    df = df[inobssample]
    #df.to_csv('RESOLVE_inobssample.csv')
    print(len(df))
    if sdsscat == 'port':
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
        heii = df['Flux_HeII_4685']
        heii_err = df['Flux_HeII_4685_Err']
    elif sdsscat =='nsa' or sdsscat == 'jhu':
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
    
    #gooddata = ((h_alpha > floor) & (nii_sum > floor) & (oiii > floor) & (oi > floor) &
    #            (sii_sum > floor) & (h_beta > floor)  & (h_alpha_err > floor) & 
    #            (nii_sum_err > floor) & (oiii_err > floor) & 
    #            (oi_err > floor) & (sii_sum_err > floor) & 
    #            (h_alpha < ceil) & (nii_sum < ceil) & (oiii < ceil) & (oi < ceil) &
    #            (sii_sum < ceil) & (h_beta < ceil) & 
    #            ~np.isnan(h_alpha) & ~np.isnan(nii_sum) & ~np.isnan(oiii) & 
    #            ~np.isnan(oi) & ~np.isnan(sii_sum) & ~np.isnan(h_beta))
    
    #df = df[gooddata]
    
    if hasnr:
        snr = ((h_alpha > cut*h_alpha_err)) 
        gooddata = ((h_alpha > floor) & (h_alpha_err > floor) & 
                (h_alpha < ceil)  & 
                ~np.isnan(h_alpha))
    
    
    elif bpt1snr:
        snr = ((h_alpha > cut*h_alpha_err) 
        & (nii_sum > cut*nii_sum_err) & 
       (oiii > cut*oiii_err) & (h_beta > cut*h_beta_err))
    
        gooddata = ((h_alpha > floor) & (nii_sum > floor) & (oiii > floor) 
                & (h_beta > floor)  & (h_alpha_err > floor) & 
                (nii_sum_err > floor) & (oiii_err > floor) & 
                (h_alpha < ceil) & (nii_sum < ceil) & (oiii < ceil) & (h_beta < ceil) 
                & ~np.isnan(h_alpha) & ~np.isnan(nii_sum) & ~np.isnan(oiii) & 
                ~np.isnan(h_beta))
    
    elif he2:
        snr = ((h_alpha > cut*h_alpha_err) 
            & (nii_sum > cut*nii_sum_err) & 
           (oiii > cut*oiii_err) & (oi > cut*oi_err) & (sii_sum > cut*sii_sum_err) 
           & (h_beta > cut*h_beta_err) &(heii > he2cut*heii_err))
        gooddata = ((h_alpha > floor) & (nii_sum > floor) & (oiii > floor) & (oi > floor) &
                (sii_sum > floor) & (h_beta > floor)  & (h_alpha_err > floor) & 
                (nii_sum_err > floor) & (oiii_err > floor) & 
                (oi_err > floor) & (sii_sum_err > floor) & (heii > floor) & (heii_err > floor) &
                (h_alpha < ceil) & (nii_sum < ceil) & (oiii < ceil) & (oi < ceil) &
                (sii_sum < ceil) & (h_beta < ceil) & 
                ~np.isnan(h_alpha) & ~np.isnan(nii_sum) & ~np.isnan(oiii) & 
                ~np.isnan(oi) & ~np.isnan(sii_sum) & ~np.isnan(h_beta) & 
                    ~np.isnan(heii))
    else:
        snr = ((h_alpha > cut*h_alpha_err) 
            & (nii_sum > cut*nii_sum_err) & 
           (oiii > cut*oiii_err) & (oi > cut*oi_err) & (sii_sum > cut*sii_sum_err) 
           & (h_beta > cut*h_beta_err))
        gooddata = ((h_alpha > floor) & (nii_sum > floor) & (oiii > floor) & (oi > floor) &
                (sii_sum > floor) & (h_beta > floor)  & (h_alpha_err > floor) & 
                (nii_sum_err > floor) & (oiii_err > floor) & 
                (oi_err > floor) & (sii_sum_err > floor) & 
                (h_alpha < ceil) & (nii_sum < ceil) & (oiii < ceil) & (oi < ceil) &
                (sii_sum < ceil) & (h_beta < ceil) & 
                ~np.isnan(h_alpha) & ~np.isnan(nii_sum) & ~np.isnan(oiii) & 
                ~np.isnan(oi) & ~np.isnan(sii_sum) & ~np.isnan(h_beta))
    print(np.sum(~np.isnan(h_alpha)))
    if resolve:
        df[~np.isnan(h_alpha)].to_csv("RESOLVE_raw_"+sdsscat+".csv")
    if eco:
        df[~np.isnan(h_alpha)].to_csv("ECO_raw_"+sdsscat+".csv")
    df = df[gooddata]
    print(len(df))
    df = df[snr]
    print(len(df)) 

    if saveoutput:
        df.to_csv(outputfile)
        print(len(df))
    
    return df

