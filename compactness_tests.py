# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 13:48:52 2021

@author: mugdhapolimera
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams.update({'axes.linewidth': 2})
matplotlib.rcParams.update({'lines.linewidth': 2})

def quenching_sfr(mass):
    return 0.64*mass-7.22

resfull = pd.read_csv('C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO+RESOLVE_barysample.csv')
resfull.index = resfull.name

resfull = resfull[resfull.logmstar > 6]
resfull['logmbary'] = np.log10(10**resfull.logmstar + 10**resfull.logmgas)

resfull = resfull[resfull.logmbary > 9.2]

resphot = pd.read_csv("Barro_inobssample.csv")
resphot.index = resphot.name 
ecophot = pd.read_csv("ECO_barysample_photometrics.csv")
ecophot.index = ecophot.name

phot = resphot.loc[np.intersect1d(resfull.name, resphot.name)]
phot = phot.append(ecophot.loc[np.intersect1d(resfull.name, ecophot.name)])

photmatch = np.intersect1d(resfull.name, phot.name)
resfull = resfull.loc[photmatch]
mass = np.arange(6.5,12.5,0.2) #np.arange(np.min(resfull.logmstar)-0.5,np.max(resfull.logmstar)+0.5,0.5)

plt.figure()
plt.plot(resfull.logmstar,np.log10(resfull.sfr_nuv),'o')
plt.plot(mass,quenching_sfr(mass))

sfg = resfull[np.log10(resfull.sfr_nuv) > quenching_sfr(resfull.logmstar)]
qg = resfull[np.log10(resfull.sfr_nuv) <= quenching_sfr(resfull.logmstar)]

plt.figure()
plt.plot(resfull.loc[photmatch].logmstar, phot.loc[photmatch].MU_50,'o')
plt.plot(resfull.loc[photmatch].loc[qg.name].logmstar,
         phot.loc[photmatch].loc[qg.name].MU_50,'ro')

qginobs = np.intersect1d(qg.name, photmatch)
compactfit = np.polyfit(np.array(resfull.loc[qginobs].logmstar), 
           np.array(phot.loc[qginobs].MU_50), 1)
plt.plot(mass,np.poly1d(compactfit)(mass),'r')
plt.plot(mass,np.poly1d(compactfit)(mass)-0.2,'k')

cg = phot[phot.MU_50 > np.poly1d(compactfit)(np.log10(phot.MSTARS))-0.2]
eg = phot[phot.MU_50 < np.poly1d(compactfit)(np.log10(phot.MSTARS))-0.2]
#cg = phot[phot.MU_DELTA > 8.6]
#eg = phot[phot.MU_DELTA < 8.6]
#cg = phot[phot.DEL_MU50 > 8.6]
#eg = phot[phot.DEL_MU50 < 8.6]

csfg = np.intersect1d(cg.name, sfg.name)
esfg = np.intersect1d(eg.name, sfg.name)

#fig, (ax1, ax2) = plt.subplots(1,2)
fig, (ax1) = plt.subplots(1,1)

#Optical
agn = pd.read_csv("ECO+RESOLVE_AGN_list.csv")
agn.index = agn.name
agn = agn[(agn.agntype != 'xrayagn') & (agn.agntype != 'midiragn')] #agn[agn.bptagn | agn.agntosf | agn.bptcomposite | agn.sfingagn]
#agn = agn[agn.agntype == 'midiragn'] #agn[agn.bptagn | agn.agntosf | agn.bptcomposite | agn.sfingagn]
csfg_agn = np.intersect1d(agn.name, csfg)
esfg_agn = np.intersect1d(agn.name, esfg)

bins = np.arange(6.5,12.5,0.2)
csfg_hist = np.histogram(resfull.loc[csfg].logmstar, bins = bins)
esfg_hist = np.histogram(resfull.loc[esfg].logmstar, bins = bins)
csfg_agn_hist = np.histogram(resfull.loc[csfg_agn].logmstar, bins = bins)
esfg_agn_hist = np.histogram(resfull.loc[esfg_agn].logmstar, bins = bins)
ax1.plot(bins[:-1]+0.1, csfg_agn_hist[0]*1.0/csfg_hist[0], 'go', ms = 15,label = 'Compact SF Galaxies')
ax1.plot(bins[:-1]+0.1, esfg_agn_hist[0]*1.0/esfg_hist[0], 'bs', ms = 15, label = 'Extended SF Galaxies')
ax1.set_ylim(0.001,1.0)
ax1.set_yscale('log')
ax1.set_xticks(bins)
ax1.set_xlim(8.5,11.5)
ax1.set_xlabel(r'$ \rm log(M_*/M_\odot)$')
ax1.set_ylabel('% of Optical AGN per mass bin')
ax1.legend(loc='upper left')
ax1.axvline(x=9.5)

#Mid-IR
#agn = pd.read_csv(r'C:\Users\mugdhapolimera\github\SDSS_spectra\mid_ir\RESOLVE_WISE_AGN.csv')
#agn.index = agn.name
#csfg_agn = np.intersect1d(agn.name, csfg)
#esfg_agn = np.intersect1d(agn.name, esfg)
#
#bins = np.arange(6.5,12.5,0.2)
#csfg_hist = np.histogram(resfull.loc[csfg].logmstar, bins = bins)
#esfg_hist = np.histogram(resfull.loc[esfg].logmstar, bins = bins)
#csfg_agn_hist = np.histogram(resfull.loc[csfg_agn].logmstar, bins = bins)
#esfg_agn_hist = np.histogram(resfull.loc[esfg_agn].logmstar, bins = bins)
#ax2.plot(bins[:-1]+0.1, csfg_agn_hist[0]*1.0/csfg_hist[0], 'gp', ms = 15,label = 'Compact SF Galaxies')
#ax2.plot(bins[:-1]+0.1, esfg_agn_hist[0]*1.0/esfg_hist[0], 'bp', ms = 15, label = 'Extended SF Galaxies')
#ax2.set_ylim(0.001,1.0)
#ax2.set_yscale('log')
#ax2.set_xticks(bins)
#ax2.set_xlim(8.5,11.5)
#ax2.set_xlabel(r'$ \rm log(M_*/M_\odot)$')
#ax2.set_ylabel('% of Mid-IR AGN per mass bin')
#ax2.legend(loc='upper left')
#ax2.axvline(x=9.5)
#
