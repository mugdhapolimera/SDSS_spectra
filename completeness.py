# -*- coding: utf-8 -*-
"""
Created on Thu May 23 10:52:56 2019

@author: mugdhapolimera
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import matplotlib
matplotlib.rcParams.update({'font.size': 20})

os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')

df = pd.read_csv('C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_full_blend_dext_new.csv')
print len(df)
inobssample = (((df.grpcz >= 3000.) & (df.grpcz <= 7000.)))
#df = df[inobssample]
ineco = (130.05 < df.radeg) & (df.radeg < 237.45)
mgas = df.logmgas
mstars = df.logmstar
full_mbary = np.log10(10**mgas + 10**mstars)
print len(df)

#eco_full = df[inobssample]
eco_mbary = np.log10(10**mgas + 10**mstars)
inobssample = (((df.grpcz >= 3000.) & (df.grpcz <= 7000.)) & ineco & \
((eco_mbary > 9.2)))# | (df.absrmag < -17.33)))
eco = df[inobssample]
eco_mbary = np.log10(10**eco.logmgas + 10**eco.logmstar)

plt.figure()
plt.plot(df.absrmag, full_mbary, 'b.')
plt.plot(df.absrmag[inobssample], full_mbary[inobssample], 'r.')
plt.plot(eco.absrmag, eco_mbary, 'ko', mfc ='none')
yaxis = np.linspace(8,12)
plt.plot(-17.33*np.ones(len(yaxis)), yaxis, 'r-.') 
xaxis = np.linspace(-24,-17)
plt.plot(xaxis, 9.2*np.ones(len(xaxis)), 'r-.')
plt.xlim(-16.5,-24)
plt.ylim(8,12)
n_eco = float(len(eco))
n_full = len(full_mbary)

print n_eco, n_full, n_eco/n_full
lowlum = np.sum((eco.absrmag > -17.33) & (eco_mbary > 9.2))
ecoabovem = float(np.sum((eco_mbary > 9.2)))
print lowlum, ecoabovem, lowlum/ecoabovem


inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_blend_dext_new.pkl'
df = pd.read_pickle(inputfile)
print len(df)
ra=df.radeg
dec=df.dedeg
grpcz = df.grpcz
cz = df.cz
infall = (ra > 22*15.) | (ra < 3*15.)
inspring = (ra > 8.75*15.) & (ra < 15.75*15.)
inobssample = ((grpcz >= 4500.) & (grpcz <= 7000.) & infall) 
df = df[inobssample]
flinsample = df.fl_insample
mgas = df.logmgas
mstars = df.logmstar
mbary = 10**mgas + 10**mstars
resfullmbary = np.log10(mbary)
inobssample = (resfullmbary > 9.2)
#inobssample = (((flinsample | (np.log10(mbary) > 9.2)) & infall))
#               | ((flinsample | (np.log10(mbary) > 9.2)) & inspring))
#inobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & \
#(((flinsample | (np.log10(mbary) > 9.0)) & infall) | \
#        ((flinsample | (np.log10(mbary) > 9.2)) & inspring))
resolve = df[inobssample]
mgas = resolve.logmgas
mstars = resolve.logmstar
mbary = 10**mgas + 10**mstars
resmbary = np.log10(mbary)

res_el = pd.read_csv('RESOLVE_full_hasnr5_dext_jhu.csv')
resmbary_el = np.log10(10**res_el.logmstar + 10**res_el.logmgas)
res_sel = pd.read_csv('RESOLVE_full_snr5_dext_jhu.csv')
resmbary_sel = np.log10(10**res_sel.logmstar + 10**res_sel.logmgas)
massthresh = 9.2
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(df.absrmag, resfullmbary, '.', color = 'gray', alpha = 0.2, 
         label = 'Parent Sample')
ax1.plot(resolve.absrmag, resmbary, '.', color = 'gray',
#         mfc ='none', 
         label = 'Mass-Limited Sample')
ax1.plot(res_el.absrmag, resmbary_el, '.', color = 'orange',
#         mfc ='none', 
         label = 'Emission Line Sample')
ax1.plot(res_sel.absrmag, resmbary_sel, '.', color = 'k',
#         mfc ='none', 
         label = 'Strong Emission Line Sample')
yaxis = np.linspace(-24,-15)
ax1.plot(yaxis, massthresh*np.ones(len(yaxis)), 'b-.')#, 
#         label = 'Mass Limit for RESOLVE-B')
ax1.text(-21.8, 9.3, 'Mass floor', color = 'b')# fontsize = 14, color = 'k')
yaxis = np.linspace(8,12)
ax1.plot(-17.0*np.ones(len(yaxis)), yaxis, 'r-.')#, 
#         label = 'Luminosity limit for RESOLVE-B ')
#ax1.plot(-17.33*np.ones(len(yaxis)), yaxis, 'g-.', 
#         label = 'Luminosity limit for ECO ')
ax1.text(-16.5, 11.5, 'Luminosity floor', color = 'r', rotation = 90)# fontsize = 14, color = 'k')
ax1.set_xlim(-15,-24)
ax1.set_ylim(8,12)
ax1.set_xlabel(r'$\rm M_r$')#, fontsize = 15)
ax1.set_ylabel(r'$\rm log(M_{baryonic}/M_\odot)$')#, fontsize = 15)
ax1.legend(fontsize = 15, loc = 4)
ax1.text(-21.8, 11.75, 'RESOLVE - B', color = 'k')# fontsize = 14, color = 'k')

lowlum = (resolve.absrmag > -17.0) & (resmbary > massthresh)
highlum = (resolve.absrmag < -17.0) & (resmbary > massthresh)
resbelowm = resmbary < massthresh
resabovem = float(np.sum((resmbary > massthresh)))
#ax1.plot(resolve.absrmag[lowlum], resmbary[lowlum], 'c.')
#ax1.plot(resolve.absrmag[highlum], resmbary[highlum], 'm.')
##ax1.plot(resolve.absrmag[resbelowm & (resolve.absrmag < -17.0)], 
##                         resmbary[resbelowm & (resolve.absrmag < -17.0)], 'g.')
#ax1.plot(resolve.absrmag[resbelowm], resmbary[resbelowm], 'g.')

massthresheco = 9.35
lowlumeco = (resolve.absrmag > -17.33) & (resmbary > massthresheco)
resaboveecom = float(np.sum((resmbary > massthresheco)))

print "RESOLVE Incompleteness"
print np.sum(lowlum), resabovem, 100.0*np.sum(lowlum)/resabovem
print "ECO Incompleteness"
print np.sum(lowlumeco), resaboveecom, 100.0*np.sum(lowlumeco)/resaboveecom
#print np.sum(lowlum), len(resolve), 100.0*np.sum(lowlum)/len(resolve)
#print np.sum(lowlum), len(resfullmbary), 100.0*np.sum(lowlum)/len(resfullmbary)

inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_snr5.pkl'
df = pd.read_pickle(inputfile)
print len(df)
mgas = df.logmgas
mstars = df.logmstar
mbary = 10**mgas + 10**mstars
finalmbary = np.log10(mbary)
resolvefinal = df
#ax2 = fig.add_subplot(122)
#ax2.plot(resolve.absrmag, resmbary, '.', color = 'k',  
#         label = 'Mass+Luminosity Limited Sample')
#ax2.plot(resolvefinal.absrmag, finalmbary, 'o', color = 'r',
#         mfc ='none', label = 'Final Sample')
##yaxis = np.linspace(-24,-15)
##ax2.plot(yaxis, 9.2*np.ones(len(yaxis)), 'b-.', 
##         label = 'Luminosity limit for ECO')
##yaxis = np.linspace(8,12)
##ax2.plot(-17.33*np.ones(len(yaxis)), yaxis, 'r-.', 
##         label = 'Mass limit for RESOLVE and ECO')
#ax2.set_xlim(-15,-24)
#ax2.set_ylim(8.5,12)
#ax2.set_xlabel(r'$\rm M_r$', fontsize = 14)
#ax2.set_ylabel(r'$\rm M_{baryonic}$', fontsize = 14)
#ax2.legend(fontsize = 14, loc = 2)
#ax2.text(-22, 8.75, 'RESOLVE - B', fontsize = 14, color = 'k')
lowlum = np.sum((resolve.absrmag > -17.33) & (resmbary > 9.2))
resabovem = float(np.sum((resmbary > 9.2)))
print lowlum, resabovem, 100.0*lowlum/resabovem
