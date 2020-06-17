# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 13:54:14 2020

@author: mugdhapolimera
Compare ECO and RESOLVE databases
"""

import pandas as pd
import numpy as np

res = pd.read_csv('RESOLVE_snr5_master.csv')
res.index = res.name
resflag = pd.read_csv('resolve_emlineclass_full_snr5_master.csv')
resflag.index = resflag.galname
resfull = pd.read_csv('RESOLVE_full_blend_dext_new.csv')
resfull.index = resfull.name

eco = pd.read_csv('ECO_snr5_master.csv')
eco.index = eco.name
ecoflag = pd.read_csv('eco_emlineclass_full_snr5_master.csv')
ecoflag.index = ecoflag.galname
ecofull = pd.read_csv('ECO_full_blend_dext_new.csv')
ecofull.index = ecofull.name

nsaeco = pd.read_csv('NSA_ECO.csv')
nsaeco.index = nsaeco.resname
nsaecofull = ecofull.loc[list(nsaeco.resname)]
lines = list(nsaeco.keys())
lines.remove('Name')
lines.remove('Name.1')
lines.remove('dedeg')
lines.remove('radeg')
lines.remove('resname')
print(lines)
nsaecofull[lines] = nsaeco[lines]
nsaecofull.to_csv('NSA_ECO_full.csv')
#extra = [x for x in list(eco.resname) if (x not in list(res.name)) and (x != 'notinresolve')]
#
#xndx = [x for x in range(len(eco)) if eco.iloc[x].resname in extra]
#df = ecofull.loc[list(eco.iloc[xndx].name)]
#nii = df['nii_6584_flux']
#if 'nii_6548_flux' in df.keys():
#    nii_sum = (df['nii_6584_flux']+ df['nii_6548_flux'])*3./4
#    nii_sum_err = (np.sqrt(df['nii_6584_flux_err']**2 + df['nii_6548_flux_err']**2))*3./4
#else:
#    nii_sum = df['nii_6584_flux']
#    nii_sum_err = df['nii_6584_flux_err']
## note that the ratio uses only the stronger line, but for S/N reasons we add
## the weaker and multiply by 3/4 since Chris Richardson says the canonical
## line ratio is 3:1 (this needs to be updated with a more precise number)
#oiii = df['oiii_5007_flux']
#oiii_err = df['oiii_5007_flux_err']
#h_alpha = df['h_alpha_flux']
#h_alpha_err = df['h_alpha_flux_err']
#h_beta = df['h_beta_flux']
#h_beta_err = df['h_beta_flux_err']
#oi = df['oi_6300_flux']
#oi_err = df['oi_6300_flux_err']
#if 'sii_6717_flux' in df.keys():
#    sii_sum = df['sii_6717_flux'] + df['sii_6731_flux']
#
#    sii_sum_err = np.sqrt(df['sii_6717_flux_err']**2 + df['sii_6731_flux_err']**2)
#else:
#    sii_sum = df['sii_6731_flux']
#
#    sii_sum_err = df['sii_6731_flux_err']
#floor = 10**-3
#ceil = 1e5
#
#gooddata = ((h_alpha > floor) & (nii_sum > floor) & (oiii > floor) & (oi > floor) &
#            (sii_sum > floor) & (h_beta > floor)  & (h_alpha_err > floor) & 
#            (nii_sum_err > floor) & (oiii_err > floor) & 
#            (oi_err > floor) & (sii_sum_err > floor) & 
#            (h_alpha < ceil) & (nii_sum < ceil) & (oiii < ceil) & (oi < ceil) &
#            (sii_sum < ceil) & (h_beta < ceil) & 
#            ~np.isnan(h_alpha) & ~np.isnan(nii_sum) & ~np.isnan(oiii) & 
#            ~np.isnan(oi) & ~np.isnan(sii_sum) & ~np.isnan(h_beta))
#
#print(len(gooddata), np.sum(gooddata))
#
