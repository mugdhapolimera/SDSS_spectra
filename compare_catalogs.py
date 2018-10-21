# -*- coding: utf-8 -*-
"""
Created on Mon Oct 01 11:59:16 2018

@author: mugdhapolimera
"""

inputfile = 'RESOLVE_SDSS_dext.fits'
dat = Table.read(inputfile, format='fits')
df1 = dat.to_pandas()

nii_sum1 = (df1['nii_6584_flux_ext'] + df1['nii_6548_flux_ext'])*3./4
oiii1 = df1['oiii_5007_flux_ext']
h_alpha1 = df1['h_alpha_flux_ext']
h_beta1 = df1['h_beta_flux_ext']
oi1 = df1['oi_6300_flux_ext']
sii_sum1 = df1['sii_6717_flux_ext'] + df1['sii_6731_flux_ext']
heii1 = df1['Flux_HeII_4685_ext']
#need to check Kewley paper and figure out if ratio used in nii_sum applies to sii_sum as well

#filter out bad data and merge gooddata with blue E/S0's
gooddata1 = ((h_alpha1 > 0) & (nii_sum1 > 0) & (oiii1 > 0) & (oi1 > 0) &
            (sii_sum1 > 0) & (h_beta1 > 0) ) & (heii > 0))
count= 0
for i in range(len(df)):
    for j in range(len(df1)):
        if df1['NAME'][j] == df['NAME'][i]:
            if abs(heii[i] - heii1[j]) >= 10**-2:
                count += 1                
                print df1['NAME'][j]           
                print df['h_beta_flux_ext'][i],df1['h_beta_flux_ext'][j]
                print df['h_alpha_flux_ext'][i],df1['h_alpha_flux_ext'][j]
                print df['oiii_5007_flux_ext'][i],df1['oiii_5007_flux_ext'][j]
                print df['nii_6584_flux_ext'][i],df1['nii_6584_flux_ext'][j]
                print df['nii_6548_flux_ext'][i], df1['nii_6548_flux_ext'][j]
                print df['oi_6300_flux_ext'][i],df1['oi_6300_flux_ext'][j]
                print df['sii_6717_flux_ext'][i], df1['sii_6717_flux_ext'][j] 
                print df['sii_6731_flux_ext'][i], df1['sii_6731_flux_ext'][j]
print count
