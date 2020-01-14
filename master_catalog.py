# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 18:59:10 2019

@author: mugdhapolimera

Making a master catalog from MPA-JHU, Portsmouth, NSA catalogs
"""

import pandas as pd
import os
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.io.idl import readsav
import numpy as np

os.chdir('C:\Users\mugdhapolimera\github\SDSS_spectra/')
resolve = pd.read_csv('RESOLVE_full_raw.csv', index_col = 'name')

jhuflag = pd.read_csv('resolve_emlineclass_full_snr5_jhu.csv', index_col = 'galname')
portflag = pd.read_csv('resolve_emlineclass_full_snr5_port.csv', index_col = 'galname')
nsaflag = pd.read_csv('resolve_emlineclass_full_snr5_nsa.csv', index_col = 'galname')

port = pd.read_csv('RESOLVE_full_snr5_port.csv', index_col = 'name')#[portflag.sftoagn]
jhu = pd.read_csv('RESOLVE_full_snr5.csv', index_col = 'name')#[jhuflag.sftoagn]
nsa = pd.read_csv('NSA_RESOLVE.csv', index_col = 'resname')#.loc[nsaflag.index.values[nsaflag.sftoagn]]
 

allunq = np.unique(list(jhuflag.index) + list(nsaflag.index) + list(portflag.index))   
#sftoagn = df[ambigsel1 & dwarf][['radeg','dedeg']]

#sfagn = pd.read_csv('uniquesfagn.csv')
unique = np.unique(list(jhuflag.index[jhuflag.sftoagn]) + \
                   list(nsaflag.index[nsaflag.sftoagn]) + \
                   list(port.index[portflag.sftoagn]))

#print('List of Unique SFing-AGN from JHU, Portsmouth and NSA Catalogs')
#print(len(unique))
#print (jhu.loc[[x for x in jhu.index if x in unique]])
#print (port.loc[[x for x in port.index \
#                 if (x in unique) & (x not in jhu.index)]])
#print (nsa.loc[[x for x in nsa.index \
#                if (x in unique) & (x not in jhu.index) & (x not in port.index)]])

df = jhu
df['source'] = 'jhu'
nsakeys = ['nii_6584_flux', 'nii_6584_flux_err', 
           'sii_6731_flux', 'sii_6731_flux_err',
           'oi_6300_flux', 'oi_6300_flux_err', 
           'oiii_5007_flux', 'oiii_5007_flux_err', 
           'h_beta_flux', 'h_beta_flux_err',
           'h_alpha_flux', 'h_alpha_flux_err']
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

jhuagn = jhuflag.sftoagn | jhuflag.composite | jhuflag.defagn
nsaagnflag = nsaflag.sftoagn | nsaflag.composite | nsaflag.defagn
portagnflag = portflag.sftoagn | portflag.composite | portflag.defagn
nsaagn = nsaflag.index.values[nsaagnflag]
portagn = portflag.index.values[portagnflag]

for gal in allunq:
    flag = True
    if gal in jhuagn.index.values:
        flag = (~jhuagn.loc[gal]) | (gal in portagn) | (gal in nsaagn)
    if flag:
        if ((gal in port.index.values) & (gal not in nsaagn)) | \
        ((gal in portagn) & (gal in nsaagn)): 
            if (gal not in df.index.values):
                print 'port', gal
                df = df.append(resolve.loc[gal])
            df['source'][gal] = 'port'
            for key in range(len(jhukeys)):
                df.set_value(gal, jhukeys[key], port.loc[gal][portkeys[key]])
        elif gal in nsaflag.index.values:
            if (gal not in df.index.values):
                print 'nsa', gal
                df = df.append(resolve.loc[gal])
            df['source'][gal] = 'nsa'
            df.loc[gal, nsakeys] = nsa.loc[gal][nsakeys]
            df.loc[gal,'sii_6717_flux'] = 0.000001
            df.loc[gal,'sii_6717_flux_err'] = 0.000001
            df.loc[gal,'nii_6548_flux'] = 0.000001
            df.loc[gal,'nii_6548_flux_err'] = 0.000001
#df.to_csv('RESOLVE_snr5_master.csv')
master = df