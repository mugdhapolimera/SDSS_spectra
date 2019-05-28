# -*- coding: utf-8 -*-
"""
Created on Thu May 23 10:52:56 2019

@author: mugdhapolimera
"""

import pandas as pd
import matplotlib.pyplot as ax2
import os
import numpy as np

os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')

df = pd.read_pickle('ECO_full_blend_dext_new.pkl')
print len(df)
inobssample = (((df.grpcz >= 3000.) & (df.grpcz <= 7000.)))
df = df[inobssample]
mgas = df.logmgas
mstars = df.logmstar
full_mbary = np.log10(10**mgas + 10**mstars)
print len(df)

#eco_full = df[inobssample]
eco_mbary = np.log10(10**mgas + 10**mstars)
inobssample = (eco_mbary > 9.2)
eco = df[inobssample]
eco_mbary = np.log10(10**eco.logmgas + 10**eco.logmstar)

#ax2.figure()
#ax2.plot(df.absrmag, full_mbary, 'b.')
#ax2.plot(eco.absrmag, eco_mbary, 'ko', mfc ='none')
#yaxis = np.linspace(-24,-17)
#ax2.plot(yaxis, 9.2*np.ones(len(yaxis)), 'r-.')
#ax2.xlim(-17,-24)
#ax2.ylim(8,12)
#n_eco = float(len(eco))
#n_full = len(full_mbary)
#
#print n_eco, n_full, n_eco/n_full

inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_blend_dext_new.pkl'
df = pd.read_pickle(inputfile)
print len(df)
ra=df.radeg
dec=df.dedeg
grpcz = df.grpcz
cz = df.cz
infall = (ra > 22*15.) | (ra < 3*15.)
inobssample = ((grpcz >= 4500.) & (grpcz <= 7000.) & infall) 
df = df[inobssample]

flinsample = df.fl_insample
mgas = df.logmgas
mstars = df.logmstar
mbary = 10**mgas + 10**mstars
resfullmbary = np.log10(mbary)
inobssample = (resfullmbary > 9.2)
#inobssample = (((flinsample | (np.log10(mbary) > 9.0)) & infall) | 
#        ((flinsample | (np.log10(mbary) > 9.2)) & inspring))
resolve = df[inobssample]
mgas = resolve.logmgas
mstars = resolve.logmstar
mbary = 10**mgas + 10**mstars
resmbary = np.log10(mbary)

fig = plt.figure()
ax1 = fig.add_subplot(121)
ax1.plot(df.absrmag, resfullmbary, '.', color = 'gray', alpha = 0.3, 
         label = 'Parent Sample')
ax1.plot(resolve.absrmag, resmbary, 'k.', 
         label = 'Mass Limited Sample')
yaxis = np.linspace(-24,-15)
ax1.plot(yaxis, 9.2*np.ones(len(yaxis)), 'b-.', 
         label = 'Luminosity limit for ECO')
yaxis = np.linspace(8,12)
ax1.plot(-17.33*np.ones(len(yaxis)), yaxis, 'r-.', 
         label = 'Mass limit for RESOLVE and ECO')
ax1.set_xlim(-15,-24)
ax1.set_ylim(8,12)
ax1.set_xlabel(r'$\rm M_r$', fontsize = 14)
ax1.set_ylabel(r'$\rm M_{baryonic}$', fontsize = 14)
ax1.legend(fontsize = 14, loc = 2)
ax1.text(-22, 8.25, 'RESOLVE - B', fontsize = 14, color = 'k')
lowlum = np.sum((resolve.absrmag > -17.33) & (resmbary > 9.2))
resabovem = float(np.sum((resmbary > 9.2)))
print lowlum, resabovem, lowlum/resabovem

inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_filter_new.pkl'
df = pd.read_pickle(inputfile)
print len(df)
mgas = df.logmgas
mstars = df.logmstar
mbary = 10**mgas + 10**mstars
finalmbary = np.log10(mbary)
resolvefinal = df
ax2 = fig.add_subplot(122)
ax2.plot(resolve.absrmag, resmbary, '.', color = 'k',  
         label = 'Mass Limited Sample')
ax2.plot(resolvefinal.absrmag, finalmbary, 'o', color = 'r',
         mfc ='none', label = 'Final Sample')
#yaxis = np.linspace(-24,-15)
#ax2.plot(yaxis, 9.2*np.ones(len(yaxis)), 'b-.', 
#         label = 'Luminosity limit for ECO')
#yaxis = np.linspace(8,12)
#ax2.plot(-17.33*np.ones(len(yaxis)), yaxis, 'r-.', 
#         label = 'Mass limit for RESOLVE and ECO')
ax2.set_xlim(-15,-24)
ax2.set_ylim(9,12)
ax2.set_xlabel(r'$\rm M_r$', fontsize = 14)
ax2.set_ylabel(r'$\rm M_{baryonic}$', fontsize = 14)
ax2.legend(fontsize = 14, loc = 2)
ax2.text(-22, 9.25, 'RESOLVE - B', fontsize = 14, color = 'k')
lowlum = np.sum((resolve.absrmag > -17.33) & (resmbary > 9.2))
resabovem = float(np.sum((resmbary > 9.2)))
print lowlum, resabovem, lowlum/resabovem
