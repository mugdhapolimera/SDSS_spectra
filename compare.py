# -*- coding: utf-8 -*-
"""
Created on Sat Sep 29 01:20:38 2018

@author: mugdhapolimera
"""
import pandas as pd
import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np

inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_SDSS_all_dext.fits'
dat = Table.read(inputfile, format='fits')
df = dat.to_pandas()
inputfile = 'C:/Users/mugdhapolimera/github/BPT/BPT/BPT/RESOLVE_SDSS_dext.fits'
dat = Table.read(inputfile, format='fits')
df1 = dat.to_pandas()

nii_sum = (df['nii_6584_flux_ext'] + df['nii_6548_flux_ext'])*3./4
oiii = df['oiii_5007_flux_ext']
h_alpha = df['h_alpha_flux_ext']
h_beta = df['h_beta_flux_ext']
h_beta_err = df['h_beta_flux_ext_err']
oi = df['oi_6300_flux_ext']
sii_sum = df['sii_6717_flux_ext'] + df['sii_6731_flux_ext']
heii = df['heii_4685_flux_port_ext']
heii_err = df['heii_4685_flux_port_ext_err']

nii_sum1 = (df1['nii_6584_flux_ext'] + df1['nii_6548_flux_ext'])*3./4
oiii1 = df1['oiii_5007_flux_ext']
h_alpha1 = df1['h_alpha_flux_ext']
h_beta1 = df1['h_beta_flux_ext']
oi1 = df1['oi_6300_flux_ext']
sii_sum1 = df1['sii_6717_flux_ext'] + df1['sii_6731_flux_ext']
heii1 = df1['Flux_HeII_4685_ext']
gooddata = ((h_alpha > 0) & (nii_sum > 0) & (oiii > 0) & (oi > 0) &
            (sii_sum > 0) & (h_beta > 0) & (heii > 0))# & (h_beta > 3*h_beta_err)) # & (heii > 3*heii_err))

gooddata1 = ((h_alpha1 > 0) & (nii_sum1 > 0) & (oiii1 > 0) & (oi1 > 0) & 
    (sii_sum1 > 0) & (h_beta1 > 0) & (heii1 > 0) )#& (h_beta > 3*h_beta_err))

count = 0
for name in df1['NAME'][gooddata1]:
    if name not in list(df['NAME'][gooddata]):
        count += 1
        print count,name
        #print df.loc[np.where(df['NAME'] == name)[0]]

'''
plt.figure(1)
plt.title('OIII')

plt.figure(2)
plt.title('H_alpha')

plt.figure(3)
plt.title('H_beta')

plt.figure(4)
plt.title('OI')

plt.figure(5)
plt.title('SII_Sum')

plt.figure(6)
plt.title('HeII')

plt.figure(7)
plt.title('NII Sum')
count = 0
for i in range(len(df)):
    for j in range(len(df1)):
        if df1['NAME'][j] == df['NAME'][i]:
            count+=0            
            plt.figure(7)
            plt.plot(nii_sum[i]-nii_sum1[j],'bo')
            plt.figure(1)
            plt.plot(oiii[i]-oiii1[j],'bo')
            plt.figure(2)
            plt.plot(h_alpha[i]-h_alpha1[j],'bo')
            plt.figure(3)
            plt.plot(h_beta[i]-h_beta1[j],'bo')
            plt.figure(4)
            plt.plot(oi[i]-oi1[j],'bo')
            plt.figure(5)                        
            plt.plot(sii_sum[i]-sii_sum1[j],'bo')
            plt.figure(6)
            plt.plot(heii[i]-heii1[j],'bo')            
'''