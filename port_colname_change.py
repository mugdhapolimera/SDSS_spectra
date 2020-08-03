# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 17:29:56 2020

@author: mugdhapolimera
"""

import pandas as pd

portkeys = ['Flux_NII_6547', 'Flux_NII_6547_Err',
            'Flux_NII_6583', 'Flux_NII_6583_Err', 
           'Flux_SII_6716', 'Flux_SII_6716_Err',
           'Flux_SII_6730', 'Flux_SII_6730_Err', 
           'Flux_OI_6300', 'Flux_OI_6300_Err', 
           'Flux_OIII_5006', 'Flux_OIII_5006_Err',
           'Flux_Hb_4861', 'Flux_Hb_4861_Err', 
           'Flux_Ha_6562', 'Flux_Ha_6562_Err']
jhukeys = ['nii_6548_flux', 'nii_6548_flux_err',
           'nii_6584_flux', 'nii_6584_flux_err', 
           'sii_6717_flux', 'sii_6717_flux_err',
           'sii_6731_flux', 'sii_6731_flux_err',
           'oi_6300_flux', 'oi_6300_flux_err', 
           'oiii_5007_flux', 'oiii_5007_flux_err', 
           'h_beta_flux', 'h_beta_flux_err',
           'h_alpha_flux', 'h_alpha_flux_err']

jhu = pd.read_csv('ECO_full_SDSS_raw.csv')
jhu.index = jhu.name
port = jhu.copy()

for gal in list(jhu.name):
    for key in range(len(jhukeys)):
        port.set_value(gal, jhukeys[key], jhu.loc[gal][portkeys[key]])

port.to_csv('ECO_full_raw_port.csv')