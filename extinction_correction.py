# -*- coding: utf-8 -*-
"""
Created on Sun Dec 02 18:13:57 2018

@author: mugdhapolimera
"""

##deredden SDSS flux measurements for RESOLVE
##October 7, 2016

##############################

import numpy as np
import pandas as pd
import extinction

import os
#pylab.ion()

#flag = 'smc'
flag = 'mw'
#os.chdir('C:\Users\mugdhapolimera\github\SDSS_Spectra')
#os.chdir('./github/SDSS_Spectra')
#open the file
data = pd.read_pickle('ECO_full_raw.pkl')
wavelengths = {'oii_3726_flux' : 3726, 'oii_3729_flux' : 3729, 
               'neiii_3869_flux' : 3869, 'h_delta_flux' : 4101, 
               'h_gamma_flux' : 4340, 'oiii_4363_flux' : 4363, 
               'h_beta_flux' : 4861, 'oiii_4959_flux': 4959, 
               'oiii_5007_flux' : 5007, 'hei_5876_flux' : 5876, 
               'oi_6300_flux' : 6300, 'nii_6548_flux' : 6548, 
               'h_alpha_flux' : 6564, 'nii_6584_flux' : 6584, 
               'sii_6717_flux' : 6717, 'sii_6731_flux' : 6731,
               'ariii7135_flux' : 7135}

''', 'Flux_HeII_3203' : 3203, 
               'Flux_OII_3726' : 3726, 'Flux_OII_3728' : 3728, 
               'Flux_NeIII_3868' : 3868, 'Flux_Hd_4101' : 4101, 
               'Flux_Hg_4340' : 4340, 'Flux_OIII_4363' : 4363, 
               'Flux_HeII_4685' : 4685, 'Flux_ArIV_4711' : 4711, 
               'Flux_Hb_4861' : 4861, 'Flux_OIII_4958' : 4958, 
               'Flux_OIII_5006' : 5006, 'Flux_HeI_5875' : 5875, 
               'Flux_OI_6300' : 6300, 'Flux_NII_6547' : 6547,
               'Flux_Ha_6562' : 6562, 'Flux_NII_6583' : 6583,
               'Flux_SII_6716' : 6716, 'Flux_SII_6730' : 6731}
'''
if flag =='mw':
    extinction_func = extinction.mwextinction
    outputfile = 'ECO_full_mwdext_new.pkl'
if flag =='smc':
    extinction_func = extinction.smcextinction
    outputfile = 'ECO_full_smcdext_new.pkl'
Alambda_Ebv= {}    
for line in wavelengths.keys():               
    Alambda_Ebv[line] = extinction_func(1/(wavelengths[line]*10**-4))
              
########reddening correct###############
#intrinsic balmer decrement
k_Hbeta = extinction_func(1/(wavelengths['h_beta_flux']*10**-4))
k_Halpha = extinction_func(1/(wavelengths['h_alpha_flux']*10**-4))
#intrinsic balmer decrement
balm_dec_int = 2.86 
#observed balmer decrement
balm_dec_obs = data['h_alpha_flux']/data['h_beta_flux']
balm_dec_obs_err=balm_dec_obs*((data['h_alpha_flux_err']/data['h_alpha_flux'])**2 + 
                    (data['h_beta_flux_err']/data['h_beta_flux'])**2)**0.5

#relation from balmer decrement to color excess
EBV_excess= (2.5/(k_Hbeta - k_Halpha))*(np.log10(balm_dec_obs/balm_dec_int)) 
EBV_excess_err= (2.5/(k_Hbeta - k_Halpha)) * (0.434/balm_dec_int) * \
                (balm_dec_obs_err/balm_dec_obs) 

###################REST FRAME###################
data_dext = data.copy()

#find conversion array for each galaxy at each wavelength
mw_A_atmwlam  = np.zeros(len(EBV_excess))
#lines = [line for line in data.keys() if 'flux' in line.lower()]
#if 'flux_ratio_2c' in lines:
#    lines.pop(np.where('flux_ratio_2c' in lines)[0])
#print lines
for galname in data.index.values:
    print galname
    #factor to relate each mw lam and color excess    
    #ext['Al_Ebv']*EBV_excess[i] #for particular galaxy
    #there is a mw_A_atmwlam for each galaxy, this index is the same for all
    
    #if np.isfinite(EBV_excess.loc[galname]): 
        #Extinction for each line in a galaxy is given by
        # A_lambda_gal = A_lambda/E(B-V) * E(B-V)_gal
    Alambda_gal = {line : Alambda_Ebv[line]*EBV_excess.loc[galname] \
                            for line in wavelengths.keys()}
    Alambda_gal_err = {line : Alambda_Ebv[line]*EBV_excess_err.loc[galname] \
                            for line in wavelengths.keys()}
            #output is dextincted avg flux for each line in each galaxy
    for line in wavelengths.keys():
         
         data_dext[line].loc[galname] = (data[line].loc[galname]*
                                     10.0**(Alambda_gal[line]/2.5))
         X = 10**(Alambda_gal[line]/2.5)
         X_err = 2.303 * X * Alambda_gal_err[line]
#         if line[0] == 'F':
#             line_err = data[line+'_Err'].loc[galname]/data[line].loc[galname]
#             data_dext[line+'_Err'].loc[galname] = data_dext[line].loc[galname]*(line_err**2 + (X_err/X)**2)**0.5
 #        else:
         line_err = data[line+'_err'].loc[galname]/data[line].loc[galname]
         data_dext[line+'_err'].loc[galname] = data_dext[line].loc[galname]*(line_err**2 + (X_err/X)**2)**0.5
                                 
#print data_dext.loc[galname]         
data_dext.to_pickle(outputfile)
#####################################################################
