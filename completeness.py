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
from scipy import stats
matplotlib.rcParams.update({'font.size': 25})
matplotlib.rcParams.update({'axes.linewidth': 2})
matplotlib.rcParams.update({'lines.linewidth': 2})

os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')

df = pd.read_csv('C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_full_blend_dext_new.csv')
print len(df)
czinobssample = (((df.grpcz >= 3000.) & (df.grpcz <= 7000.)))
#df = df[inobssample]
ineco = (((130.05 < df.radeg) & (df.radeg < 237.45)) & ((-1 < df.dedeg) & (df.dedeg < 49.85)))
mgas = df.logmgas
mstars = df.logmstar
full_mbary = np.log10(10**mgas + 10**mstars)
print len(df)
#eco_full = df[inobssample]
eco_mbary = np.log10(10**mgas + 10**mstars)
inobssample = (czinobssample & ineco & ((eco_mbary > 9.2)))# | (df.absrmag < -17.33)))
eco = df[inobssample]
eco_mbary = np.log10(10**eco.logmgas + 10**eco.logmstar)
eco.to_csv("ECO_barysample.csv")

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
inobssample = ((grpcz >= 4500.) & (grpcz <= 7000.) & (infall | inspring)) 
df = df[inobssample]
parent = df
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
resolve.to_csv("RESOLVE_barysample.csv")
mgas = resolve.logmgas
mstars = resolve.logmstar
mbary = 10**mgas + 10**mstars
resmbary = np.log10(mbary)

full = resolve.copy()
notinresolve = (eco['resname'] == 'notinresolve')
full = full.append(eco[notinresolve])
full.index = full.name


res_el = pd.read_csv('ECO+RESOLVE_hasnr5_dext_jhu.csv')
resmbary_el = np.log10(10**res_el.logmstar + 10**res_el.logmgas)
res_sel = pd.read_csv('ECO+RESOLVE_snr5_dext_jhu.csv')
resmbary_sel = np.log10(10**res_sel.logmstar + 10**res_sel.logmgas)
massthresh = 9.2
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(df.absrmag, resfullmbary, '.', color = 'gray', alpha = 0.2, 
         label = 'Full Sample')
ax1.plot(resolve.absrmag, resmbary, '.', color = 'black',
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
ax1.text(-16.5, 11.75, 'RESOLVE-B Luminosity floor', color = 'r', rotation = 90)# fontsize = 14, color = 'k')
ax1.plot(-17.33*np.ones(len(yaxis)), yaxis, 'g-.') 
ax1.text(-17.5, 11.75, 'ECO Luminosity floor', color = 'g', rotation = 90)# fontsize = 14, color = 'k')
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

massthresheco = 9.2
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


#bins = np.arange(-24,-15, 0.5)
#plt.figure()
#plt.hist(parent.absrmag, normed = True, bins = bins, color = 'gainsboro', 
#         label = 'Parent Sample')
#plt.hist(resolve.absrmag, normed = True, histtype = 'step', hatch = 'X', bins = bins,\
#         color = 'dimgray', lw = 2, label = 'Mass-limited Sample')
#plt.hist(res_el.absrmag, normed = True, histtype = 'step', bins = bins,
#         color = 'orange', lw = 2, label = 'Emission Line Sample')
#plt.hist(res_sel.absrmag, normed = True, histtype = 'step', bins = bins,
#         color = 'black', lw = 2, label = 'Strong Emission Line Sample')
#plt.axvline(-17.0)
#plt.xlim(-15,-24)
#plt.xlabel(r"$\rm M_r$", fontsize = 15)
#plt.ylabel("Normalized # of galaxies")
#plt.legend(loc = "upper right", fontsize = 15)
#
#bins = np.arange(7.0,12.0, 0.2)
#plt.figure()
#plt.hist(parent.logmstar, normed = True, bins = bins, color = 'gainsboro', 
#         label = 'Parent Sample')
#plt.hist(resolve.logmstar, normed = True, histtype = 'step', hatch = 'X', bins = bins,\
#         color = 'dimgray', lw = 2, label = 'Mass-limited Sample')
#plt.hist(res_el.logmstar, normed = True, histtype = 'step', bins = bins,
#         color = 'orange', lw = 2, label = 'Emission Line Sample')
#plt.hist(res_sel.logmstar, normed = True, histtype = 'step', bins = bins,
#         color = 'black', lw = 2, label = 'Strong Emission Line Sample')
##plt.axvline(9.2)
#plt.xlim(6.5,12)
#plt.xlabel(r"$\rm M_{stellar}$", fontsize = 15)
#plt.ylabel("Normalized # of galaxies")
#plt.legend(loc = "upper right", fontsize = 15)

resolve = full.copy()
resmbary = np.log10(10**resolve.logmstar + 10**resolve.logmgas)
resnsa = pd.read_csv('C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO+RESOLVE_snr5_dext_nsa.csv')
resport = pd.read_csv('C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO+RESOLVE_snr5_dext_port.csv')
resnsa_mbary = np.log10(10**resnsa.logmstar + 10**resnsa.logmgas)
resport_mbary = np.log10(10**resport.logmstar + 10**resport.logmgas)
bins = np.arange(9.2,11,0.2)
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)
#plt.hist(np.log10(mbary), normed = True, bins = bins, color = 'gainsboro', 
#         label = 'Parent Sample')
portline = ax1.hist(resport_mbary, normed = True, histtype = 'step', bins = bins,
         color = 'red', ls = '--',lw = 3, label = 'Portsmouth SEL Sample',
         zorder = 5)
nsaline = ax1.hist(resnsa_mbary, normed = True, histtype = 'step', bins = bins,
         color = 'blue', lw = 3, ls = '-.',label = 'NSA SEL Sample',
         zorder = 6)
ax1.legend(fontsize = 15)
ax1.hist(resmbary, normed = True, bins = bins,\
         color = 'lightblue', alpha = 0.8, lw = 2, label = 'Mass-limited Parent Sample')
ax1.hist(resmbary_el, normed = True, histtype = 'step', bins = bins,
         color = 'orange', lw = 4, label = 'EL Sample')
ax1.hist(resmbary_sel, normed = True, histtype = 'step', bins = bins,
         color = 'black', lw = 3, hatch = 'X', label = 'MPA-JHU SEL Sample')
#plt.axvline(9.2)
ax1.set_xlim(9.0,11.2)
ax1.set_xlabel(r"$\rm log(M_{baryonic}/M_\odot)$")
ax1.set_ylabel("Normalized # of galaxies")
#plt.legend(loc = "upper right", fontsize = 15)
ax1.set_xticks(bins)
ax1.tick_params(labelsize = 14)
ax2.tick_params(labelsize = 14)
ax3.tick_params(labelsize = 14)

bins = np.arange(10.5,14,0.25)
#plt.hist(parent.logmh, normed = True, bins = bins, color = 'gainsboro', 
#         label = 'Parent Sample')
ax2.hist(resolve.logmh, normed = True, bins = bins,\
          color = 'lightblue', alpha = 0.8, lw = 2, label = 'Mass-limited Parent Sample')
ax2.hist(res_el.logmh, normed = True, histtype = 'step', bins = bins,
         color = 'orange', lw = 4, label = 'Emission Line Sample')
ax2.hist(res_sel.logmh, normed = True, histtype = 'step', bins = bins,
         color = 'black', lw = 3, hatch = 'X', label = 'Strong Emission Line Sample')
#plt.axvline(9.2)
ax2.set_xlim(10.25,14.)
ax2.set_xticks(np.arange(10.5,14,0.5))
ax2.set_xlabel(r"$\rm log(M_{halo}/M_\odot)$")
#ax2.set_ylabel("Normalized # of galaxies")
#ax2.legend(loc = "upper right", fontsize = 15)
ax1.tick_params(length=8, width=2)
ax2.tick_params(length=8, width=2)

bins = np.arange(0.5,3.5,0.25)
ax3.hist(resolve.modelu_rcorr, normed = True, bins = bins,\
          color = 'lightblue', alpha = 0.8, lw = 2, label = 'Mass-limited Parent Sample')
ax3.hist(res_el.modelu_rcorr, normed = True, histtype = 'step', bins = bins,
         color = 'orange', lw = 4, label = 'MPA-JHU EL Sample')
ax3.hist(res_sel.modelu_rcorr, normed = True, histtype = 'step', bins = bins,
         color = 'black', lw = 3, hatch = 'X', label = 'MPA-JHU SEL Sample')
#plt.axvline(9.2)
ax3.set_xlim(0.25,3.5)
ax3.set_xticks(np.arange(0.5,3.25,0.5))
ax3.set_xlabel(r"$\rm (u-r)^e$ color")
#ax3.set_ylabel("Normalized # of galaxies")
ax3.legend(loc = "upper right", fontsize = 15)
ax3.tick_params(length=8, width=2)

resolve_oi = np.log10(resolve.oi_6300_flux)#/resolve.h_alpha_flux)
res_el_oi = np.log10(res_el.oi_6300_flux)#/res_el.h_alpha_flux)
res_sel_oi = np.log10(res_sel.oi_6300_flux)#/res_sel.h_alpha_flux)
plt.figure()
#bins = np.arange(-2.0,0.2,0.2)
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

bins = np.arange(-2.0,5.0,0.5)
ax1.hist(resolve_oi, normed = False, bins = bins,\
          color = 'lightblue', alpha = 0.8, lw = 2, label = 'Mass-limited Parent Sample')
ax1.hist(res_el_oi, normed = False, histtype = 'step', bins = bins,
         color = 'orange', lw = 4, label = 'MPA-JHU EL Sample')
ax1.hist(res_sel_oi, normed = False, histtype = 'step', bins = bins,
         color = 'black', lw = 3, hatch = 'X', label = 'MPA-JHU SEL Sample')
#plt.axvline(9.2)
ax1.set_xlim(-2.1,4.6)
ax1.set_xticks(np.arange(-2.0,4.5,1.0))
#plt.xlabel(r"log([O I])/H$ \rm \alpha$")
ax1.set_xlabel("log([O I])")
ax1.set_ylabel("# of galaxies")
#ax1.legend(loc = "upper right", fontsize = 15)
ax1.tick_params(length=8, width=2)

resolve_oiha = np.log10(resolve.oi_6300_flux/resolve.h_alpha_flux)
res_el_oiha = np.log10(res_el.oi_6300_flux/res_el.h_alpha_flux)
res_sel_oiha = np.log10(res_sel.oi_6300_flux/res_sel.h_alpha_flux)
#plt.figure()
bins = np.arange(-2.0,0.2,0.2)
#bins = np.arange(-2.0,5.0,0.5)
ax2.hist(resolve_oiha, normed = False, bins = bins,\
          color = 'lightblue', alpha = 0.8, lw = 2, label = 'Mass-limited Parent Sample')
ax2.hist(res_el_oiha, normed = False, histtype = 'step', bins = bins,
         color = 'orange', lw = 4, label = 'MPA-JHU EL Sample')
ax2.hist(res_sel_oiha, normed = False, histtype = 'step', bins = bins,
         color = 'black', lw = 3, hatch = 'X', label = 'MPA-JHU SEL Sample')
#plt.axvline(9.2)
ax2.set_xlim(-2.1,0.3)
ax2.set_xticks(np.arange(-2.0,0.2,0.5))
ax2.set_xlabel(r"log([O I]/H$ \rm \alpha)$")
#plt.xlabel("log([O I])")
#ax2.set_ylabel("# of galaxies")
ax2.legend(loc = "upper right", fontsize = 17)
ax2.tick_params(length=8, width=2)

plt.figure()
plt.scatter(resolve_oi, resolve_oiha, color = 'lightblue', lw = 2, label = 'Mass-limited Parent Sample')
plt.scatter(res_el_oi, res_el_oiha, color = 'orange', lw = 4, label = 'MPA-JHU EL Sample')
plt.scatter(res_sel_oi, res_sel_oiha, color = 'black', lw = 3, label = 'MPA-JHU SEL Sample')

#DD, pnullks = stats.ks_2samp(resmbary,resmbary_el)
#print('K-S test p-value = ' + str(pnullks))
#
#DD, pnullks = stats.ks_2samp(resmbary, resmbary_sel)
#print('K-S test p-value = ' + str(pnullks))
#
#DD, pnullks = stats.ks_2samp(resolve.logmh, res_el.logmh)
#print('K-S test p-value = ' + str(pnullks))
#
#DD, pnullks = stats.ks_2samp(resolve.logmh, res_sel.logmh)
#print('K-S test p-value = ' + str(pnullks))
#
#import scipy.stats as st
#samp_a = np.sort(resmbary_el)
#samp_b = np.sort(resmbary_sel)
#samp_conc = np.sort(np.concatenate((samp_a, samp_b)))
#samp_a_cdf = [np.round(st.percentileofscore(samp_a, value)/100, 1) for value in samp_conc]
#
##cdf of sample b
#samp_b_cdf = [np.round(st.percentileofscore(samp_b, value)/100, 1) for value in samp_conc]
##print('samp_b_cdf:\n{}'.format(samp_b_cdf))
##print(20*'-')
#
##compute absolute difference
#samp_diff = np.abs(np.subtract(samp_a_cdf, samp_b_cdf))
##print('samp_diff:\n{}'.format(samp_diff))
#plt.figure(); 
#plt.plot(samp_conc,samp_a_cdf)
#plt.plot(samp_conc,samp_b_cdf)
##plt.plot(samp_conc,min(samp_a_cdf,samp_b_cdf)+samp_diff)
#
#print("D_crit ", 1.36*np.sqrt(1.0/len(samp_a) + 1.0/len(samp_b)))
#print("D_n ", max(samp_diff))
