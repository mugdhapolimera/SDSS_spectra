# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 13:30:34 2022

@author: mugdhapolimera
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 14:54:24 2019

@author: mugdhapolimera

Mid-IR AGN Selection using WISE IR colours

as prescribed by Sartori et al. 2015 using the method of Jarrett et al 2015, 
Stern et al 2012; and the method of Satyapal et al. 2018 

GAMA data from 
http://www.gama-survey.org/dr3/schema/table.php?id=56
"""

import numpy as np
import pandas as pd
pd.set_option('display.max_rows',500)
import matplotlib.pyplot as plt
import os
from scipy.io.idl import readsav
from astropy.io import fits
#os.chdir('C:\Users\mugdhapolimera\github\SDSS_spectra')
def stern(x):
    return 0.8*np.ones(len(x))

def jarretx(y):
    return [2.2*np.ones(len(y)), 4.2*np.ones(len(y))]

def jarrety(x):
    return [1.7*np.ones(len(x)), 0.1*x+0.38]

def satyapalx(x):
    return 0.52*np.ones(len(x))

def satyapaly(x):
    return 5.78*x -24.50
##############################################################################
#Reading in RESOLVE catalogs and new WISE photometry and setting up the data
##############################################################################
def midir_plotallagn(ax, midirfile, s06outputfile, seloutputfile, survey):
    df = pd.read_csv(midirfile)
    df.index = df.name
    df = df[~np.isnan(df.mw1)]
    sel = pd.read_csv(seloutputfile)
    sel.index = sel.galname
    selagn = np.array(sel[~sel.defstarform].galname)

    s06 = pd.read_csv(s06outputfile)
    s06.index = s06.galname
    s06agn = np.array(s06[~s06.defstarform].galname)

    df_full = df.copy()
    
    print('Total '+survey+' galaxies: {}'.format(len(df)))
    
    ##############################################################################
    #Performing quality control on the data
    ##############################################################################
    
    #Removing nans, 0 and applying S/N > 5 thresholding for the mid-IR mags
    baderr = np.isnan(df.emw1) | np.isnan(df.emw2) | np.isnan(df.emw3)
    badphot = (df.mw1 <= 0.0) | (df.mw2 <= 0.0) | (df.mw3 <= 0.0) 
    threshold = 5
    fluxw1=10**(df.mw1/(-2.5))
    efluxw1=df.emw1*(np.log10(10.0)/(2.5))*fluxw1
    fluxw2=10**(df.mw2/(-2.5))
    efluxw2=df.emw2*(np.log10(10.0)/(2.5))*fluxw2
    fluxw3=10**(df.mw3/(-2.5))
    efluxw3=df.emw3*(np.log10(10.0)/(2.5))*fluxw3
    fluxw4=10**(df.mw4/(-2.5))
    efluxw4=df.emw4*(np.log10(10.0)/(2.5))*fluxw4
    goodwise1 = ((df.mw1 > 0) & (fluxw1/efluxw1 > 5.0))
    goodwise2 = ((df.mw2 > 0) & (fluxw2/efluxw2 > 5.0))
    goodwise3 = ((df.mw3 > 0) & (fluxw3/efluxw3 > 5.0)) 
    snr = goodwise1 & goodwise1 & goodwise3

    good = ~baderr & ~badphot & snr
    print('Galaxies with true mags/errors and S/N > {}: {}'.format(threshold, \
          np.sum(good)))
        
    ##############################################################################
    #AGN classification based on WISE colour magnitudes
    ##############################################################################
    w12 = df['mw1'] - df['mw2']
    w23 = df['mw2'] - df['mw3']
    w12_err = np.sqrt(df['emw1']**2 + df['emw2']**2)
    w23_err = np.sqrt(df['emw2']**2 + df['emw3']**2)
    
    midiragn = ((w12 >= 0.8) | ((w23 > 2.2) & (w23 < 4.2) & \
                 (w12< 1.7) & (0.1*w23 + 0.38 < w12)) | \
               (w12>= 0.52) & (w12>= (5.78*w23) -24.50))

    #plt.figure()
    xaxis = np.linspace(min(w23)-0.1,max(w23)+0.1)
    #yaxis = np.linspace(min(w12), max(w12))
    yaxis = np.linspace(jarrety(np.array([2.2]))[1],1.7)
    ax.plot(xaxis, stern(xaxis), 'k-.')#, label = 'Stern12')
    ax.text(-0.8,0.9,'St12', fontsize = 15)
    ax.plot(xaxis, satyapalx(xaxis), 'k')#, label = 'Satyapal18')
    ax.text(-0.8,0.6 ,'Sa18', fontsize = 15)
    xaxis = np.linspace(4.3287,max(w23))
    ax.plot(xaxis, satyapaly(xaxis), 'k')
    
    xaxis = np.linspace(2.2,4.2)
    ax.plot(jarretx(yaxis)[0], yaxis, 'k--', jarretx(yaxis)[1], yaxis, 'k--')
    ax.plot(xaxis, jarrety(xaxis)[0], 'k--')
    ax.plot(xaxis, jarrety(xaxis)[1],'k--')#, label = 'Jarrett15')
    ax.text(4,1.5,'J11', fontsize = 15)
    ax.set_xlabel('W2 - W3')
    ax.set_ylabel('W1 - W2')
    ax.set_ylim(min(w12)-0.1, max(w12)+0.1)
#    ax.plot(w23['rs0775'],w12['rs0775'],'p', color = 'black', ms = 10, mec = 'none',
#             label = 'rs0775')
    dwarfs = df.logmstar < 9.5
    dwarfagn = dwarfs & midiragn
    ax.plot(w23[midiragn],w12[midiragn],'p', color = 'orange', ms = 10, mec = 'none', zorder = 2.,
             label = 'Mid-IR AGN')
    ax.plot(w23[dwarfagn],w12[dwarfagn],'kp', ms = 12, mec = 'k', mfc = 'none',
                 label = 'Mid-IR Dwarf AGN', zorder = 12)
    ax.set_ylim(-1.2, 2.2)
    ax.set_xlim(-1.2,4.5)#min(w23)-0.1,max(w23))
    
    seldwarfagn = np.intersect1d(selagn, df.name[dwarfs])
    ax.plot(w23[seldwarfagn],w12[seldwarfagn],'ks', ms = 12, mec = 'k', mfc = 'none',
                 label = 'SEL Dwarf AGN', zorder = 10)
    s06dwarfagn = np.intersect1d(s06agn, df.name[dwarfs])
    ax.plot(w23[s06dwarfagn],w12[s06dwarfagn],'kv', ms = 12, mec = 'k', mfc = 'none',
                 label = 'S06 Dwarf AGN', zorder = 11)
#    ax.plot(w23['rf0206'],w12['rf0206'],'rs', ms = 12, mec = 'r', mfc = 'none',
#                 label = 'S06 Dwarf AGN', zorder = 11)
    ax.legend(loc = 'upper left', fontsize = 15)
    ax.plot(w23,w12,'gp', alpha = 0.3, ms = 10, mec = 'none', zorder = 0,
             label = 'Mid-IR SF')
    
    midiragnname = df.name[midiragn]
    #print(resolve.loc[midiragnname][['radeg','dedeg']])    
    df_midiragn = df.loc[midiragnname]
    print(list(midiragnname))
    df['agnflag'] = False
    df['agnflag'][midiragn] = True
    
    print('{} mid-IR AGN out of {} galaxies having reliable WISE mags : {}%'\
          .format(len(midiragnname), len(df), \
                  round(len(midiragnname)*100.0/len(df),2)))
    print('Dwarf mid-IR galaxies: ',np.sum(dwarfs))
    print('Dwarf mid-IR AGN: ',np.sum(dwarfagn))
    
