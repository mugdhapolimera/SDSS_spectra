# -*- coding: utf-8 -*-
"""
Created on Mon May 20 15:56:52 2019

@author: mugdhapolimera

RESOLVE and ECO data visualization plots
"""

import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')
df = pd.read_pickle('RESOLVE_full_blend_dext_new.pkl')
#df = pd.read_csv('RESOLVE_snr5_master.csv')
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
inobssample = ((grpcz >= 4500.) & (grpcz <= 7000.)) & (((flinsample | (np.log10(mbary) > 9.0)) & infall) | ((flinsample | (np.log10(mbary) > 9.2)) & inspring))
#inobssample = (((grpcz >= 4500.) & (grpcz <= 7000.)) & \
#((((np.log10(mbary) > 9.0)) & infall) | \
# (((np.log10(mbary) > 9.2)) & inspring)))
resolve_full = df[inobssample]

df = pd.read_pickle('ECO_full_blend_dext_new.pkl')
mgas = df.logmgas
mstars = df.logmstar
mbary = 10**mgas + 10**mstars
inobssample = (((df.grpcz >= 3000.) & (df.grpcz <= 7000.)) & \
((np.log10(mbary) > 9.2))) #(df.absrmag < -17.3) & 
eco_full = df[inobssample]


resolve = pd.read_csv('RESOLVE_snr5_master.csv')
eco = pd.read_csv('ECO_filter_new.csv')

resolve_mbary_full = np.log10(10**resolve_full.logmstar +  
                            10**resolve_full.logmgas)
eco_mbary_full = np.log10(10**eco_full.logmstar +  10**eco_full.logmgas)
resolve_mbary = np.log10(10**resolve.logmstar +  10**resolve.logmgas)
eco_mbary = np.log10(10**eco.logmstar +  10**eco.logmgas)

fig = plt.figure('Baryonic Mass Distributions')
ax1 = fig.add_subplot(111)
ax1.hist(resolve_mbary_full, bins = 'fd', ec = 'orange', 
         histtype = 'step', lw = 5, label = 'Parent RESOLVE Sample \n(1519 galaxies)')
ax1.hist(resolve_mbary, bins = 'fd', histtype = 'step', 
         hatch = 'x', ec = 'k', lw = 5 , label = 'Final Sample \n(437 galaxies)')
ax1.set_ylabel('Number of Galaxies', size = 15)
ax1.set_xlabel(r'$\rm \log (M_{baryonic})$ $\rm M_{\odot}$', size = 15)
legend = ax1.legend(fontsize = 14) #title = 'RESOLVE'
#legend.get_title().set_fontsize('14')

#ax1 = fig.add_subplot(122)
#ax1.hist(eco_mbary_full, bins = 'fd', histtype = 'step', lw = 5, 
#         label = 'Parent Sample')
#ax1.hist(eco_mbary, bins = 'fd', histtype = 'stepfilled', 
#         hatch = 'x', ec = 'k', lw = 5, label = 'Final Sample')
#ax1.set_ylabel('Number of Galaxies', size = 15)
#ax1.set_xlabel(r'$\rm \log (M_{baryonic})$ $\rm M_{\odot}$', size = 15)
#legend = ax1.legend(title = 'ECO', fontsize = 14)
#legend.get_title().set_fontsize('14')

resolve_full_coords = SkyCoord(ra=list(resolve_full.radeg)*u.degree, 
                           dec=list(resolve_full.dedeg)*u.degree, frame='icrs')
resolve_coords = SkyCoord(ra=list(resolve.radeg)*u.degree, 
                      dec=list(resolve.dedeg)*u.degree, frame='icrs')

eco_full_coords = SkyCoord(ra=list(eco_full.radeg)*u.degree, 
                           dec=list(eco_full.dedeg)*u.degree, frame='icrs')
eco_coords = SkyCoord(ra=list(eco.radeg)*u.degree, 
                      dec=list(eco.dedeg)*u.degree, frame='icrs')

fig = plt.figure('Spatial Distributions')
ax1 = fig.add_subplot(121)
ax1.scatter(list(resolve_full_coords.ra.hour), list(resolve_full.dedeg), 
            c = 'orange', alpha = 0.3, label = 'Parent Sample')
ax1.plot(resolve_coords.ra.hour, resolve.dedeg, 'k.', label = 'Final Sample')
ax1.set_ylabel('Declination (degrees)', size = 15)
ax1.set_xlabel('Right Ascension (hours)', size = 15)
legend = ax1.legend(title = 'RESOLVE', fontsize = 14, loc = 2)
legend.get_title().set_fontsize('14')
ax1.set_xlim(-5, 25)

ax2 = fig.add_subplot(122)
ax2.scatter(list(eco_full_coords.ra.hour), list(eco_full.dedeg), 
            c = 'orange', alpha = 0.3, label = 'Parent Sample')
ax2.plot(eco_coords.ra.hour, eco.dedeg, 'k.')
ax2.set_ylabel('Declination (degrees)', size = 15)
ax2.set_xlabel('Right Ascension (hours)', size = 15)
ax2.set_xlim(7.75,16)
legend = ax2.legend([], [], title = 'ECO', loc = 2)
legend.get_title().set_fontsize('14')

#resolve = pd.read_csv('ECO+RESOLVE_filter_new.csv')
#N2 = np.log10(resolve.nii_6584_flux/resolve.h_alpha_flux)
#selgd<-which(N2>(-2.5) & N2<(-0.3))
#selbd = ((N2<-2.5) | (N2>-0.3))
#resolve_Z = 9.37+2.03*N2+1.26*N2**2+0.32*N2**3
zdata = np.genfromtxt('RESOLVE+ECO_bpassagn_nichollsjen_ssp20.txt', 
        dtype = None, names = ["name", "Estimate", "err_up", "err_down"]) 
resolve_Z = zdata["Estimate"]+8.76
#rep = np.isnan(resolve_Z)
#resolve_Z[rep] = (-99.)
#resolve_Z[selbd] = (-99.)
resolve_bpt = pd.read_csv('eco+resolve_emlineclass_filter.csv')
resolve_sfagn_Z = resolve_Z[list(resolve_bpt['sftoagn'])]
resolve_defagn_Z = resolve_Z[list(resolve_bpt['defagn'])]

N2 = np.log10(eco.nii_6584_flux/eco.h_alpha_flux)
#selgd<-which(N2>(-2.5) & N2<(-0.3))
selbd = ((N2<-2.5) | (N2>-0.3))
eco_Z = 9.37+2.03*N2+1.26*N2**2+0.32*N2**3
rep = np.isnan(eco_Z)
eco_Z[rep] = (-99.)
eco_Z[selbd] = (-99.)
eco_bpt = pd.read_csv('eco_emlineclass_filter_new.csv')
eco_sfagn_Z = eco_Z[list(eco_bpt['sftoagn'])]
eco_defagn_Z = eco_Z[list(eco_bpt['defagn'])]


fig = plt.figure('Metallicity Distributions')
ax1 = fig.add_subplot(121)
weights = np.ones_like(resolve_Z)/float(len(resolve_Z))
hist_full = ax1.hist(resolve_Z[resolve_Z > -50], bins = 'fd', histtype = 'step',  
         ec = 'green', lw = 3, density = True, label = 'All Galaxies')
hist_sfagn = ax1.hist(resolve_sfagn_Z, bins = 'fd', histtype = 'step', 
            hatch = 'x', ec = 'orange', lw = 3, density = True, 
            label = 'SF-to-AGN Galaxies')
#hist_defagn = ax1.hist(resolve_defagn_Z, bins = 'fd', histtype = 'step', 
#            hatch = '\/', ec = 'blue', lw = 3, density = True, 
#            label = 'SF-to-AGN Galaxies')
ax1.set_ylabel('Number of Galaxies', size = 15)
ax1.set_xlabel(r'Metallicity ($\rm 12 + \log (O/H)$ using BPASS model)', size = 15)

histmax = max(max(hist_full[0]), max(hist_sfagn[0]))
x_solar = np.linspace(0, 1.5*histmax)
y_solar = 8.76*np.ones(len(x_solar))
y_solar_40 = 8.76+np.log10(0.4)*np.ones(len(x_solar))
y_solar_30 = 8.76+np.log10(0.3)*np.ones(len(x_solar))
y_solar_20 = 8.76+np.log10(0.2)*np.ones(len(x_solar))

ax1.plot(y_solar_40,x_solar,'k-.', linewidth = 2)
ax1.plot(y_solar_30,x_solar,'k-.', linewidth = 2)
ax1.plot(y_solar_20,x_solar,'k-.', linewidth = 2)
ax1.plot(y_solar,x_solar,'k-.', linewidth = 2)
ax1.text(y_solar_20[0]-0.05, 1.2*histmax, r'0.2 $Z_{\odot}$',
         fontsize=14, color='k', rotation = 'vertical')
ax1.text(y_solar_30[0]-0.05, 1.2*histmax, r'0.3 $Z_{\odot}$',
         fontsize=14, color='k', rotation = 'vertical')
ax1.text(y_solar_40[0]-0.05, 1.2*histmax, r'0.4 $Z_{\odot}$',
         fontsize=14, color='k', rotation = 'vertical')
ax1.text(y_solar[0]-0.05, 1.2*histmax, r'1.0 $Z_{\odot}$',
         fontsize=14, color='k', rotation = 'vertical')
legend = ax1.legend(title = 'RESOLVE + ECO', loc = 2)
legend.get_title().set_fontsize('14')
ax1.set_ylim(0,max(x_solar))

ax1 = fig.add_subplot(122)
weights = np.ones_like(eco_Z)/float(len(eco_Z))
hist_full = ax1.hist(eco_Z[eco_Z > -50], bins = 'fd', histtype = 'step',  
         ec = 'green', lw = 3, density = True, label = 'All Galaxies')
hist_sfagn = ax1.hist(eco_sfagn_Z, bins = 'fd', histtype = 'step', 
        hatch = 'x', ec = 'orange', lw = 3, density = True, 
        label = 'SF-to-AGN Galaxies')
#hist_defagn = ax1.hist(eco_defagn_Z, bins = 'fd', histtype = 'step', 
#        hatch = '\/', ec = 'blue', lw = 3, density = True, 
#        label = 'SF-to-AGN Galaxies')
ax1.set_ylabel('Number of Galaxies', size = 15)
ax1.set_xlabel(r'Metallicity (in $\rm 12 + \log (O/H)$ units)', size = 15)

histmax = max(max(hist_full[0]), max(hist_sfagn[0]))
x_solar = np.linspace(0, 1.5*histmax)
y_solar = 8.76*np.ones(len(x_solar))
y_solar_40 = 8.76+np.log10(0.4)*np.ones(len(x_solar))
y_solar_30 = 8.76+np.log10(0.3)*np.ones(len(x_solar))
y_solar_20 = 8.76+np.log10(0.2)*np.ones(len(x_solar))

ax1.plot(y_solar_40,x_solar,'k-.', linewidth = 2)
ax1.plot(y_solar_30,x_solar,'k-.', linewidth = 2)
ax1.plot(y_solar_20,x_solar,'k-.', linewidth = 2)
ax1.plot(y_solar,x_solar,'k-.', linewidth = 2)
ax1.text(y_solar_20[0]-0.05, 1.2*histmax, r'0.2 $Z_{\odot}$',
         fontsize=14, color='k', rotation = 'vertical')
ax1.text(y_solar_30[0]-0.05, 1.2*histmax, r'0.3 $Z_{\odot}$',
         fontsize=14, color='k', rotation = 'vertical')
ax1.text(y_solar_40[0]-0.05, 1.2*histmax, r'0.4 $Z_{\odot}$',
         fontsize=14, color='k', rotation = 'vertical')
ax1.text(y_solar[0]-0.05, 1.2*histmax, r'1.0 $Z_{\odot}$',
         fontsize=14, color='k', rotation = 'vertical')
legend = ax1.legend(title = 'ECO', loc = 2)
legend.get_title().set_fontsize('14')
ax1.set_ylim(0,max(x_solar))
