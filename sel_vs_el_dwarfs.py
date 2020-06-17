# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 15:59:39 2020

@author: mugdhapolimera

Compare SEL and EL dwarfs
"""
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 15:25:30 2018

@author: mugdhapolimera

This code explores the properties of galaxies categorized using BPT plots.
"""

import numpy as np
import os
import pandas as pd
pd.set_option('display.max_columns', 500)
from matplotlib.ticker import NullFormatter
from matplotlib.ticker import ScalarFormatter
#from scipy.stats import norm

#path = os.path.dirname(os.path.realpath(__file__))
os.chdir('C:/Users/mugdhapolimera/github/SDSS_Spectra/')
he2_flag = 0
full = 0
resolve = 1
eco = 0
#if he2_flag:
#    flags = pd.read_csv('eco+resolve_emlineclass_filter_he2.csv')
#else:
#    flags = pd.read_csv('eco+resolve_emlineclass_filter.csv')
if full : 
    inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO+RESOLVE_filter_new.pkl'
    flags = pd.read_csv('eco+resolve_emlineclass_filter.csv')
if resolve: 
    inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_hasnr3_master.csv'
    flags = pd.read_csv('resolve_emlineclass_full_snr5_master.csv')
if eco: 
    inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_filter_new.pkl'
    flags = pd.read_csv('eco_emlineclass_filter_new.csv')

full_df = pd.read_csv(inputfile, index_col = 'name')
full_df = full_df[full_df.logmstar < 9.5]
df = pd.read_csv('RESOLVE_snr5_master.csv', index_col = 'name')
df = df[df.logmstar<9.5]
#full_df.loc[list(flags['galname'])]

    
#G/S vs. M_star with histogram
'''rect_histy = [left_h, bottom, 0.2, height]
axScatter = plt.axes(rect_scatter)
axHisty = plt.axes(rect_histy, xscale = 'log')
nullfmt = NullFormatter() # no labels
axHisty.yaxis.set_major_formatter(nullfmt)
'''
plt.figure('Gas Content')
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left + width + 0.02
rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.2]
rect_histy = [left_h, bottom, 0.2, height]
axScatter = plt.axes(rect_scatter)
axHistx = plt.axes(rect_histx, yscale = 'log')
nullfmt = NullFormatter() # no labels
axHistx.xaxis.set_major_formatter(nullfmt)

# the scatter plot:
axScatter.plot(full_df.logmstar,
                   10**full_df.logmgas/10**full_df.logmstar, 
                   'k.', markersize = 8, 
                   alpha = 0.3, mew = 0, label = 'EL Dwarfs')

axScatter.plot(df.logmstar,
                   10**df.logmgas/10**df.logmstar, 
                   'bs', markersize = 12, mew = 0, 
                   label = 'SEL Dwarfs')
axScatter.plot(np.linspace(7.5,12.5), 
               np.ones(len(np.linspace(7.5,12.5))), 'k-.')
#axScatter.text(11.0, 0.005, 'RESOLVE', fontsize=14, color='k')
axScatter.text(10.5, 1.1, r'1:1 G/S Ratio', fontsize=14, color='k')

axScatter.text(9.52, -2.2, r'Gas Richness', fontsize=14, color='k')
axScatter.text(9.52, -2.4, r'Threshold Mass', fontsize=14, color='k')

axScatter.plot(9.5*np.ones(len(np.linspace(10**-2.5,10**1.5))), 
               np.linspace(10**-2.5,10**1.5), 'k-.')
axScatter.set_ylim(10**-2.5,10**1.5)
axScatter.set_xlim(7.5,12.5) 
axScatter.set_yscale("log")
axScatter.yaxis.set_major_formatter(ScalarFormatter())
axScatter.legend(loc='upper right',numpoints = 1, fontsize = 12) # bbox_to_anchor=(1,1),
#if he2_flag:
    #axScatter.legend(loc=2, bbox_to_anchor=(1,1.15),numpoints = 1)

axScatter.set_xlabel(r'$\rm \log(M_{stellar}/M_{\odot})$', fontsize = 22)
axScatter.set_ylabel(r'$\rm M_{gas}/M_{stellar}$',fontsize = 22)

axHistx.hist(full_df.logmstar, histtype = 'stepfilled', alpha = 0.1,
         bins= bins, linewidth = 5, label = 'EL Dwarfs')
axHistx.hist(df.logmstar, histtype = 'step',bins = bins, alpha = 0.3,
             linewidth = 5, label = 'SEL Dwarfs', 
                color = 'blue')
axHistx.legend(loc='upper right', fontsize = 12) #bbox_to_anchor=(1,1.15), 
    
axHistx.set_ylabel('Number')
axHistx.set_xlim(axScatter.get_xlim())
axHistx.yaxis.set_major_formatter(ScalarFormatter())
 
plt.figure()
n,bins,patches = plt.hist(full_df.logmstar, histtype = 'stepfilled', alpha = 0.1,
         bins= 'fd', linewidth = 5, label = 'EL Dwarfs', density = True, stacked = True)
plt.hist(df.logmstar, histtype = 'step',bins = bins, density = True, stacked = True,
             linewidth = 5, label = 'SEL Dwarfs', 
                color = 'blue')
plt.legend(loc='upper right', fontsize = 12)
plt.xlabel('Stellar Mass')

plt.figure()
fullmbary = np.log10(10**full_df.logmstar + 10**full_df.logmgas)
mbary = np.log10(10**df.logmstar + 10**df.logmgas)
n,bins,patches = plt.hist(fullmbary, histtype = 'stepfilled', alpha = 0.1,
         bins= 'fd', linewidth = 5, label = 'EL Dwarfs', density = True, stacked = True)
plt.hist(mbary, histtype = 'step',bins = bins, density = True, stacked = True,
             linewidth = 5, label = 'SEL Dwarfs', 
                color = 'blue')
plt.legend(loc='upper right', fontsize = 12)
plt.xlabel('Baryonic Mass')

plt.figure()
fullg_s = full_df.logmstar - full_df.logmgas
g_s = df.logmstar - df.logmgas
n,bins,patches = plt.hist(fullg_s, histtype = 'stepfilled', alpha = 0.1,
         bins= 'fd', linewidth = 5, label = 'EL Dwarfs', density = True, stacked = True)
plt.hist(g_s, histtype = 'step',bins = bins, density = True, stacked = True,
             linewidth = 5, label = 'SEL Dwarfs', 
                color = 'blue')
plt.legend(loc='upper right', fontsize = 12)
plt.xlabel('G/S')

#FSMGR vs SSFR
plt.figure()
fullssfr = np.log10(full_df.sfr_nuv) - full_df.logmstar
ssfr = np.log10(df.sfr_nuv) - df.logmstar
plt.plot(fullssfr, full_df.meanfsmgr,
     'k.', markersize = 10, alpha = 0.3,  mew = 0, 
     label = 'EL Dwarfs')
plt.plot(ssfr, df.meanfsmgr,
     'bs', markersize = 10, mew = 0, 
      label = 'SEL Dwarfs')

#plt.plot(0*np.linspace(0.001,100), np.linspace(0.001,100), 'k-.')
plt.plot(np.linspace(-12,-8), 1+0*np.linspace(-12,-8), 'k-.')
#plt.text(-0.1, 10**-1.5, r'1:1 G/S Ratio', fontsize=14, color='k', 
#         rotation = 'vertical')
plt.text(-11.9, 1.5, r'Stellar Mass Doubled in last Gyr', 
         fontsize=14, color='k')

plt.xlabel('SSFR (Short Term SF)', size = 22)
plt.ylabel('FSMGR (Long Term SF)', size = 22)
plt.yscale('log')
yticks = plt.yticks()[0]
plt.yticks(yticks, np.around(yticks,2))
plt.ylim(10**-3, 10**2)
plt.xlim(-12,-8.5)
#plt.legend(title = 'RESOLVE', loc = 'lower right', fontsize = 14)
plt.figure()
n,bins,patches = plt.hist(fullssfr, histtype = 'stepfilled', alpha = 0.1,
         bins= 'fd', linewidth = 5, label = 'EL Dwarfs', density = True, stacked = True)
plt.hist(ssfr, histtype = 'step',bins = bins, density = True, stacked = True,
             linewidth = 5, label = 'SEL Dwarfs', 
                color = 'blue')
plt.legend(loc='upper right', fontsize = 12)
plt.xlabel('SSFR Short Term')
plt.figure()
n,bins,patches = plt.hist(full_df.meanfsmgr, histtype = 'stepfilled', alpha = 0.1,
         bins= 'fd', linewidth = 5, label = 'EL Dwarfs', density = True, stacked = True)
plt.hist(df.meanfsmgr, histtype = 'step',bins = bins, density = True, stacked = True,
             linewidth = 5, label = 'SEL Dwarfs', 
                color = 'blue')
plt.legend(loc='upper right', fontsize = 12)
plt.xlabel('FSMGR Long Term')

#Mass-Metallicity
plt.figure()
ax1 = plt.subplot()
fullN2 = np.log10(full_df.nii_6584_flux/full_df.h_alpha_flux)
N2 = np.log10(df.nii_6584_flux/df.h_alpha_flux)
bad = ((N2<-2.5) | (N2>-0.3))
Z_pp04 = 9.37+2.03*N2+1.26*N2**2+0.32*N2**3
Z = pd.read_csv("C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_snr5_master_csf_nicholls_fstar30_agn.txt", 
                     sep = '\s+', names = ["name", "Estimate", "err_up", "err_down"])
Z.index = Z.name
Z_pp04 = Z['Estimate']+8.76
#Z_pp04[bad] = -99
fullZ = (fullN2 + 15.614)/1.754
Z = (N2 + 15.614)/1.754
ax1.plot(full_df.logmstar,fullZ,
         'k.', alpha = 0.3, label = 'EL Dwarfs')
ax1.plot(df.logmstar, Z,
 'bs', markersize = 10, label = 'SEL Dwarfs')

yaxis = np.linspace(np.min(Z)-0.1,np.max(Z)+0.1)
xaxis = np.linspace(7.5, 11.5)
plt.plot(9.5*np.ones(len(yaxis)),yaxis, 'k-.', linewidth = 3)
Z04 = np.log10(0.4)+8.76
plt.plot(xaxis, Z04*np.ones(len(xaxis)), 'k--', linewidth = 2)
plt.plot(xaxis, 8.76*np.ones(len(xaxis)), 'k--', linewidth = 2)
#plt.ylim(min(yaxis),max(yaxis))
#plt.ylim(7.8,9.2)
#plt.xlim(7.5,11.5)
plt.xlabel(r'log(M${_{stellar}}$/M${_\odot}$)', fontsize = 22)
plt.ylabel(r'Z (12 + log(O/H))', fontsize = 22)
#plt.legend(fontsize = 12)
##ax2.set_xticks([0.0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,1.0])
##xticks = np.array([7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2])
#ax2 = ax1.twinx()
#yticks = np.array([7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0,9.2])
#ax2.set_yticks(np.linspace(0,1,len(yticks)))
#ax1.set_yticks(yticks)#np.arange(7.8,9.2,0.2))
#float_formatter = lambda x: "%.2f" % x
#N2 = 1.754*yticks - 15.614
#N2_label = ["%.2f" % z for z in N2]
#ax2.set_yticklabels(N2_label)
#ax2.set_ylabel(r'[NII]/H$\alpha$', fontsize = 22)

plt.figure()
n,bins,patches = plt.hist(fullZ, histtype = 'stepfilled', alpha = 0.1,
         bins= 'fd', linewidth = 5, label = 'EL Dwarfs', density = True, stacked = True)
plt.hist(Z, histtype = 'step',bins = bins, density = True, stacked = True,
             linewidth = 5, label = 'SEL Dwarfs', 
                color = 'blue')
plt.legend(loc='upper right', fontsize = 12)
plt.xlabel('Metallicity')

plt.figure()
n,bins,patches = plt.hist(full_df.logmh, histtype = 'stepfilled', alpha = 0.1,
         bins= 'fd', linewidth = 5, label = 'EL Dwarfs', density = True, stacked = True)
plt.hist(df.logmh, histtype = 'step',bins = bins, density = True, stacked = True,
             linewidth = 5, label = 'SEL Dwarfs', 
                color = 'blue')
plt.legend(loc='upper right', fontsize = 12)
plt.xlabel('Group Halo Mass from HAM')

plt.figure()
bins = np.arange(1,20)
plt.hist(full_df.grpnassoc, histtype = 'stepfilled', alpha = 0.1,
         bins= bins, linewidth = 5, label = 'EL Dwarfs', density = True, stacked = True)
plt.hist(df.grpnassoc, histtype = 'step',bins = bins, density = True, stacked = True,
             linewidth = 5, label = 'SEL Dwarfs', 
                color = 'blue')
plt.legend(loc='upper right', fontsize = 12)
plt.xlabel('Number of Galaxies in group')
