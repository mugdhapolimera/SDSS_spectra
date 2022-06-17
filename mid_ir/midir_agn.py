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
def midiragnplot(ax, inobssamplefile, survey, save):
#    resolve = pd.read_csv('C:\Users\mugdhapolimera\github\SDSS_spectra\RESOLVE_full_raw.csv')
    inobssample = pd.read_csv(inobssamplefile)
    inobssample = inobssample.set_index('name')
    inobssample = inobssample[(10**inobssample.logmstar + 10**inobssample.logmgas) > 10**9.2]
    inobsname = inobssample.index
    path = "C:\Users\mugdhapolimera\github\SDSS_Spectra\mid_ir\/" #os.getcwd()+'\/'
    os.chdir(path)
    
    #Reading in the RESOLVE internal database files
    if survey == 'RESOLVE':
        catalog = readsav(r'C:\Users\mugdhapolimera\github\SDSS_Spectra\resolvecatalogphot_031622.dat')
        catalog2 = readsav(r'C:\Users\mugdhapolimera\github\SDSS_Spectra\resolvecatalog_031622.dat')
    
        #Create a dictionary of mid-IR magnitudes, K-band magnitudes, Surface brightness
        d = {'name' : catalog2['name'],
                      'mw1' : catalog['mw1w'],
                      'mw2' : catalog['mw2w'],
                      'mw3' : catalog['mw3w'],
                      'mw4' : catalog['mw4w'],
                      'emw1' : catalog['emw1w'],
                      'emw2' : catalog['emw2w'],
                      'emw3' : catalog['emw3w'],
                      'emw4' : catalog['emw4w'],
                      'kmag' : catalog['kmag'],
                      'ekmag' : catalog['ekmag'],
                      'ukidsskmag' : catalog['ukidsskmag'],
                      'ukidsskflag' : catalog['ukidsskflag'],
                      'mur50' : catalog['mur50'],
                      'logmstar' : np.log10(catalog2['mstars'])}
        #Changing byter order from IDL dtype format to pandas default dtype
        for x in d.keys():
            d[x] = d[x].byteswap().newbyteorder() 
    else:
        catalog = readsav(r'C:\Users\mugdhapolimera\github\SDSS_Spectra\ecofull_wise_wisemask_061322.dat')
#        catalog = readsav(r'C:\Users\mugdhapolimera\github\SDSS_Spectra\ecoSEL_wise_wisemask_031022.dat')
        catalog = catalog.resolve_wise['ap'][0]
#        catalog2 = pd.read_csv(r'C:\Users\mugdhapolimera\github\SDSS_Spectra\sfr_nuv_wisew_ecoSEL_newsnr.txt')
        ecodat = readsav(r'C:\Users\mugdhapolimera\github\SDSS_Spectra\eco_wresa_032918.dat')

        match, ecodatndx, catalogndx = np.intersect1d(ecodat['names'], 
                                    catalog['name'], return_indices = True)        
        #Create a dictionary of mid-IR magnitudes, K-band magnitudes, Surface brightness
        d = {'name' : np.array(ecodat['econames'][ecodatndx]).byteswap().newbyteorder(),
              'mw1' : np.array(catalog['mw1']),
                      'mw2' : np.array(catalog['mw2']),
                      'mw3' : np.array(catalog['mw3']),
                      'mw4' : np.array(catalog['mw4']),
                      'emw1' : np.array(catalog['emw1']),
                      'emw2' : np.array(catalog['emw2']),
                      'emw3' : np.array(catalog['emw3']),
                      'emw4' : np.array(catalog['emw4']),
#                      'kmag' : catalog['kmag'],
#                      'ekmag' : catalog['ekmag'],
#                      'ukidsskmag' : catalog['ukidsskmag'],
#                      'ukidsskflag' : catalog['ukidsskflag'],
#                      'mur50' : catalog['mur50'],
                      'logmstar' : np.array(ecodat['rpgoodmstarsnew'][ecodatndx]).byteswap().newbyteorder()}

    #Making a pandas dataframe using the dictionary
    df = pd.DataFrame(data = d)
    df =  df.set_index('name')
    
#    #Reading in the new WISE photometry values and converting into a pandas DF
#    wise = readsav(path+'resolve_wise_102919.dat')
#    wisedf = pd.DataFrame.from_records(data = wise['resolve_wise'][0][0])
#    wisedf = wisedf.apply(lambda x: x.astype(str(x.dtype).replace('>','')))
#    wisedf.columns = [x.lower() for x in wisedf.columns.values]
#    wisedf.index = wisedf.name
#    
#    #Update the original DF with the new photometry
#    #SKIP if you want to use the original photometry
#    df.update(wisedf)
#    gama = pd.read_csv('GAMA_WISE_RESOLVE.csv')
#    gama.index = gama.resname
#    reliable = (gama['PHOTFLAG_W1'] > 0) & (gama['PHOTFLAG_W2'] > 0) & \
#                (gama['PHOTFLAG_W3'] > 0)
#    
#    overlap = df.loc[gama.resname]
#    gama_snr = ((gama['mw2']/gama['emw2'] > overlap['mw2']/overlap['emw2']) \
#                & (gama['mw1']/gama['emw1'] > overlap['mw1']/overlap['emw1']) \
#                & (gama['mw3']/gama['emw3'] > overlap['mw3']/overlap['emw3']))
#    res_zero = (df['mw1']==0.0) | (df['mw2']==0.0) | (df['mw3']==0.0)
#    gama_full = gama.copy()
#    gama = gama[reliable & res_zero]
    
    #df.update(gama)

    df = df.loc[inobsname]    
    df = df.merge(inobssample, on='name')
    print df.keys()
    df['logmstar'] = df['logmstar_x']
    df_full = df.copy()
    
    print('Total '+survey+' galaxies: {}'.format(len(df)))
    
    ##############################################################################
    #Performing quality control on the data
    ##############################################################################
    
    #Removing nans, 0 and applying S/N > 5 thresholding for the mid-IR mags
    baderr = np.isnan(df.emw1) | np.isnan(df.emw2) | np.isnan(df.emw3) #| \
            #np.isnan(df.emw4)
    badphot = (df.mw1 <= 0.0) | (df.mw2 <= 0.0) | (df.mw3 <= 0.0) #& \
                #(df.mw4 == 0.0)
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
    snr = goodwise1 & goodwise2 & goodwise3


#    snr = (df.mw1/df.emw1 > threshold) & (df.mw2/df.emw2 > threshold) & \
#            (df.mw3/df.emw3 > threshold) #& (df.mw4/df.emw4 > threshold)
    good = ~baderr & ~badphot & snr
    #df = df[good] #Removing bad data from the DF
    print('Galaxies with true mags/errors and S/N > {}: {}'.format(threshold, \
          np.sum(good)))
    
#    #Checking UKIDSS and 2MASS k-band magnitudes- 
#    #both k-band mags > 0; 2MASS s/n > 5
#    #df.ukidsskflag = [int(x) for x in df.ukidsskflag if x != '']
#    if 'kmag' in df.keys():
#        good_kmag = ((df['ukidsskflag'] == '   0') & (df['ukidsskmag'] > 0.0)) | \
#                (df['kmag']/df['ekmag'] > 5.0) 
#    
#    #UKIDSS and 2MASS mags should be within 10 percent of each other 
#    #since they have similar wavelength bands
#        kmag_agree = (df['ukidsskmag']/df['kmag'] > 0.90) & \
#                (df['ukidsskmag']/df['kmag'] < 1.1)
#    
#        kmagflags = good_kmag #& kmag_agree
#    #df = df[kmagflags] #Removing data with bad kmags from DF
#        print('Galaxies with reliable k-band mags: {}'.format(np.sum(kmagflags)))
#    
#    #Surface brightness cut
#    if 'mur50' in df.keys():
#        sb_threshold = 23
#        sb = df['mur50'] < sb_threshold
#    #df = df[sb]
#        print('Galaxies with Surface Brightness < {}mags/arcsec^2: {}'\
#              .format(sb_threshold,np.sum(sb)))
#        #df = df[good]# & kmagflags]# & sb]
#    fulldf = df.copy()
    df= df[good]
    
    ##############################################################################
    #AGN classification based on WISE colour magnitudes
    ##############################################################################
    w12 = df['mw1'] - df['mw2']
    w23 = df['mw2'] - df['mw3']
    w12_err = np.sqrt(df['emw1']**2 + df['emw2']**2)
    w23_err = np.sqrt(df['emw2']**2 + df['emw3']**2)
    
    #mid-IR AGN if data+/-error satisfies the AGN criteria
    #Stern OR Jarrett
    #midiragn = ((w12 >= 0.8) | ((w23 > 2.2) & (w23 < 4.2) & (w12 < 1.7) & 
    #                            (0.1*w23 + 0.38 < w12)) | 
    #           (w12 >= 0.52) & (w12 >= (5.78*w23) -24.50))
    
    #Stern OR Jarrett OR Satyapal
#    midiragn = ((w12-w12_err >= 0.8) | ((w23 > 2.2) & (w23 < 4.2) & \
#                 (w12-w12_err < 1.7) & (0.1*w23 + 0.38 < w12-w12_err)) | \
#               (w12-w12_err >= 0.52) & (w12-w12_err >= (5.78*w23) -24.50))
    #Without error consideration
    midiragn = ((w12 >= 0.8) | ((w23 > 2.2) & (w23 < 4.2) & \
                 (w12< 1.7) & (0.1*w23 + 0.38 < w12)) | \
               (w12>= 0.52) & (w12>= (5.78*w23) -24.50))
    #Stern AND Jarret AND Satyapal
    #midiragn = (((w12-w12_err) >= 0.8) & (((w23-w23_err) > 2.2) & ((w23+w23_err) < 4.2) & \
    #             (w12-w12_err < 1.7) & (0.1*(w23-w23_err) + 0.38 < w12-w12_err)) & \
    #           ((w12-w12_err) >= 0.52) & ((w12-w12_err) >= (5.78*w23) -24.50))
    
    #plt.plot(w23[midiragn],w12[midiragn],'rs')
    #plt.ylim(-6.0,10)
    
    #plt.figure()
    xaxis = np.linspace(min(w23)-0.1,max(w23)+0.1)
    #yaxis = np.linspace(min(w12), max(w12))
    yaxis = np.linspace(jarrety(np.array([2.2]))[1],1.7)
    ax.plot(xaxis, stern(xaxis), 'k-.')#, label = 'Stern12')
    ax.text(5.75,1.0,'St12', fontsize = 15)
    ax.plot(xaxis, satyapalx(xaxis), 'k')#, label = 'Satyapal18')
    ax.text(5.75,0.3,'Sa14', fontsize = 15)
    ax.text(4.7,2.0 ,'Sa18', fontsize = 15)
    xaxis = np.linspace(4.3287,max(w23))
    ax.plot(xaxis, satyapaly(xaxis), 'k')
    
    xaxis = np.linspace(2.2,4.2)
    ax.plot(jarretx(yaxis)[0], yaxis, 'k--', jarretx(yaxis)[1], yaxis, 'k--')
    ax.plot(xaxis, jarrety(xaxis)[0], 'k--')
    ax.plot(xaxis, jarrety(xaxis)[1],'k--')#, label = 'Jarrett15')
    ax.text(3.5,1.85,'J11', fontsize = 15)
    ax.set_xlabel('W2 - W3')
    ax.set_ylabel('W1 - W2')
    ax.set_ylim(min(w12)-0.1, max(w12)+0.1)
    #plt.errorbar(w23,w12,fmt = 'bo', xerr = w23_err,
    #             yerr = w12_err, label = 'Galaxies with reliale WISE mags')
    #plt.errorbar(w23[midiragn],w12[midiragn],fmt = 'rs', xerr = w23_err[midiragn],
    #             yerr = w12_err[midiragn], label = 'Mid-IR AGN')
    #plt.errorbar(w23['rs0107'],w12['rs0107'],fmt = 'ks', xerr = w23_err['rs0107'],
    #             yerr = w12_err['rs0107'], label = 'rs0107')
    #plt.plot((fulldf[~good]['mw2'] - fulldf[~good]['mw3']),
    #         (fulldf[~good]['mw1'] - fulldf[~good]['mw2']),
    #         'k.', alpha = 0.3, mec = 'none',
    #         label = 'Galaxies without reliable WISE mags')
    ax.plot(w23,w12,'gp', alpha = 0.3, ms = 10, mec = 'none',
             label = 'Mid-IR SF')
    ax.plot(w23[midiragn],w12[midiragn],'p', color = 'orange', ms = 10, mec = 'none',
             label = 'Mid-IR AGN')
#    ax.plot(w23['rs0775'],w12['rs0775'],'p', color = 'black', ms = 10, mec = 'none',
#             label = 'rs0775')
    dwarfs = df.logmstar < 9.5
    dwarfagn = dwarfs & midiragn
    ax.plot(w23[dwarfagn],w12[dwarfagn],'kp', ms = 12, mec = 'k', mfc = 'none',
                 label = 'Mid-IR Dwarf AGN')
    ax.set_ylim(-1.2, 2.2)
    ax.set_xlim(-1.25,6.5)#min(w23)-0.1,max(w23))
    
    ax.legend(loc = 'lower right', fontsize = 18)
    
    midiragnname = df.index[midiragn]
    #print(resolve.loc[midiragnname][['radeg','dedeg']])    
    df_midiragn = inobssample.loc[midiragnname]
    print(list(midiragnname))
    df['agnflag'] = False
    df['agnflag'][midiragn] = True
    
    if save:
        df_midiragn.to_csv(survey+'_WISE_AGN.csv')
        df.to_csv(survey+'_WISE_good.csv')
    print('{} mid-IR AGN out of {} galaxies having reliable WISE mags : {}%'\
          .format(len(midiragnname), len(df), \
                  round(len(midiragnname)*100.0/len(df),2)))
    print('Dwarf mid-IR galaxies: ',np.sum(dwarfs))
    print('Dwarf mid-IR AGN: ',np.sum(dwarfagn))
    
##Flags to check SFing-AGN    
#flags = pd.read_csv('../resolve_emlineclass_dext_snr5_jhu.csv')
#flags.index = flags.galname
#sfagn = list(flags.galname.iloc[np.where(flags.sftoagn)])
##sfagn = unique #list(unqsfagn)
#sfagnndx = [x for x in range(len(df.name)) if df.name[x] in sfagn]
##plt.plot(w23[sfagnndx],w12[sfagnndx],'gs', label = 'SFing-AGN')
##plt.errorbar(w23[sfagnndx],w12[sfagnndx],fmt = 'gs', xerr = w23_err[sfagnndx],
##             yerr = w12_err[sfagnndx], label = 'X-ray follow up?')
#
#sfagnmidir = [df.name[midiragn][x] for x in range(len(df.name[midiragn])) \
#              if df.name[midiragn][x] in sfagn]
#print sfagnmidir
#cat = ['jhu','port','nsa']
#for c in cat:
#    ressel = pd.read_csv("../RESOLVE_full_snr5_dext_"+c+".csv")
#    midiragn_sel = [x for x in list(df_midiragn.name) if x in list(ressel.name)]
#    df_midiragn.loc[midiragn_sel].to_csv("RESOLVE_snr5_midiragn_"+c+".csv")
#    if c == 'jhu':
#        jhuiragn = midiragn_sel
#    if c == 'port':
#        portiragn = midiragn_sel
#    if c == 'nsa':
#        nsairagn = midiragn_sel
#print("JHU")
#print jhu.loc[jhuiragn][['radeg','dedeg']]
#print("Portsmouth")
#print port.loc[portiragn][['radeg','dedeg']]
#print("NSA")
#print nsa.loc[nsairagn][['radeg','dedeg']]
#
#print("JHU")
#print jhuflag[jhu.logmstar<9.5].loc[jhuiragn]
#print("Port")
#print portflag[port.logmstar<9.5].loc[portiragn]
#print("NSA")
#print nsaflag[nsa.logmstar<9.5].loc[nsairagn]
#plt.errorbar(w23.loc[jhuiragn],w12.loc[jhuiragn],fmt = 'gs', ms= 10,xerr = w23_err.loc[jhuiragn],
#             yerr = w12_err.loc[jhuiragn], label = 'JHU Mid-IR AGN')
#plt.errorbar(w23.loc[portiragn],w12.loc[portiragn],fmt = 'ks', ms = 10,xerr = w23_err.loc[portiragn],
#             yerr = w12_err.loc[portiragn], label = 'Port Mid-IR AGN')
#plt.errorbar(w23.loc[nsairagn],w12.loc[nsairagn],fmt = 'ms', ms = 10,xerr = w23_err.loc[nsairagn],
#             yerr = w12_err.loc[nsairagn], label = 'NSA Mid-IR AGN')
#plt.errorbar(w23.loc['rs0775'],w12.loc['rs0775'],fmt = 'bs', ms = 8,xerr = w23_err.loc['rs0775'],
#             yerr = w12_err.loc['rs0775'], label = 'NSA Mid-IR AGN')
#
#plt.legend()
#
#keys = ['defstarform', 'defagn', 'composite', 'agntosf', 'sftoagn']
#    
#marker = {'agntosf': 'g^', 'composite': 'ms', 'defagn': 'rs', 
#          'defstarform': 'k.','sftoagn': 'bs'}
#
#colors = {'agntosf': 'g', 'composite': 'm', 'defagn': 'r', 
#          'defliner': 'y', 'defseyf': 'c', 'heiisel': 'k',
#          'defstarform': 'gray', 'sftoagn': 'b', 'sftoagn2': 'b'}
#
#labels = {'agntosf': 'AGN-to-SF', 'ambigagn': 'Ambiguous AGN', 
#          'composite': 'Composite', 'defagn': 'Definite AGN', 
#          'defliner': 'LINER', 'defseyf': 'Seyfert', 
#          'heiisel': 'HeII-Selected AGN', 'defstarform': 'Definite SF', 
#          'sftoagn': 'SFing-AGN', 'sftoagn2' : 'MP-AGN2'}
#rescat = pd.read_csv("C:\Users\mugdhapolimera\github\SDSS_spectra\RESOLVE_full_snr5_dext_jhu.csv")
#rescat.index = rescat.name
#plt.figure()
#xaxis = np.linspace(min(w23),max(w23))
##yaxis = np.linspace(min(w12), max(w12))
#yaxis = np.linspace(jarrety(np.array([2.2]))[1],1.7)
#plt.plot(xaxis, stern(xaxis), 'k-.', label = 'Stern12')
#xaxis = np.arange(min(w23),4.34,0.01)
#plt.plot(xaxis, satyapalx(xaxis), 'k', label = 'Satyapal18')
#xaxis = np.linspace(4.33,max(w23))
#plt.plot(xaxis, satyapaly(xaxis), 'k')
#
#xaxis = np.linspace(2.2,4.2)
#plt.plot(jarretx(yaxis)[0], yaxis, 'k--')
#yaxis = np.linspace(0.8,1.7)
#plt.plot(jarretx(yaxis)[1], yaxis, 'k--')
#plt.plot(xaxis, jarrety(xaxis)[0], 'k--')
#plt.plot(xaxis, jarrety(xaxis)[1],'k--', label = 'Jarrett15')
#plt.xlabel('W2 - W3')
#plt.ylabel('W1 - W2')
#plt.ylim(-4, 4)
#plt.xlim(2.2,4.2)
#for key in keys:
#    sel = list(rescat.name.loc[(rescat.logmstar > 9.5) & (flags[key])])
#    if key == 'defstarform':
#        plt.errorbar(w23.loc[sel],w12.loc[sel], 
#             fmt = marker[key], markersize = 10, alpha = 0.3,  mew = 0, 
#             color = colors[key], label = labels[key])
#    elif key == 'agntosf': 
#        plt.errorbar(w23.loc[sel],w12.loc[sel], 
#             fmt = marker[key], markersize = 10, mew = 1, color = colors[key],
#             mec = 'y',alpha = 0.3, label = labels[key])
#    else:
#        plt.errorbar(w23.loc[sel],w12.loc[sel],
#             fmt = marker[key], markersize = 10, mew = 0, color = colors[key],
#             alpha = 0.3, label = labels[key])
#plt.legend()
#for key in keys:
#    sel = list(rescat.name.loc[(rescat.logmstar < 9.5) & (flags[key])])
#
#    if key == 'defstarform':
#        plt.errorbar(w23.loc[sel],w12.loc[sel], 
#                     fmt = marker[key], markersize = 10, alpha = 0.9,  mew = 0, 
#             color = colors[key])
#    elif key == 'agntosf': 
#        plt.errorbar(w23.loc[sel],w12.loc[sel],
#                     fmt = marker[key], markersize = 10, mew = 1, color = colors[key],
#             mec = 'y')
#    else:
#        plt.errorbar(w23.loc[sel],w12.loc[sel],
#             fmt = marker[key], markersize = 10, mew = 0, color = colors[key])
#
#
#
#df_full['kmag_flag'] = kmagflags
#df_full['sb_flag'] = sb
#df_full['wise_good'] = good
#df_full['radeg'] = resolve['radeg']
#df_full['dedeg'] = resolve['dedeg']
##df_full.to_csv('resolve_irphot_new.csv')
#
##from astropy.table import Table as table
##gama = table.read('WISECat.fits')
##mags = ['mw1','mw2','mw3']
##for x in mags:
##    plt.figure()
##    plt.errorbar(gama[x],df_full.loc[gama.resname][x],fmt='o',
##                 xerr = gama['e'+x],yerr = df_full.loc[gama.resname]['e'+x])
##    plt.xlabel('GAMA '+x)
##    plt.ylabel('RESOLVE '+x)
##    plt.xlim(min(gama[x]),max(gama[x]))
##    plt.ylim(min(df_full.loc[gama.resname][x]),max(df_full.loc[gama.resname][x]))
#
#    d = {'name' : catalog2['name'],
#                  'mw1' : catalog['mw1o'],
#                  'mw2' : catalog['mw2o'],
#                  'mw3' : catalog['mw3o'],
#                  'mw4' : catalog['mw4o'],
#                  'emw1' : catalog['emw1o'],
#                  'emw2' : catalog['emw2o'],
#                  'emw3' : catalog['emw3o'],
#                  'emw4' : catalog['emw4o'],
#                  'kmag' : catalog['kmag'],
#                  'ekmag' : catalog['ekmag'],
#                  'ukidsskmag' : catalog['ukidsskmag'],
#                  'ukidsskflag' : catalog['ukidsskflag'],
#                  'mur50' : catalog['mur50']}
#    d = {'name' : catalog2['name'],
#                  'mw1' : catalog['mw1'],
#                  'mw2' : catalog['mw2'],
#                  'mw3' : catalog['mw3'],
#                  'mw4' : catalog['mw4'],
#                  'emw1' : catalog['emw1'],
#                  'emw2' : catalog['emw2'],
#                  'emw3' : catalog['emw3'],
#                  'emw4' : catalog['emw4'],
#                  'kmag' : catalog['kmag'],
#                  'ekmag' : catalog['ekmag'],
#                  'ukidsskmag' : catalog['ukidsskmag'],
#                  'ukidsskflag' : catalog['ukidsskflag'],
#                  'mur50' : catalog['mur50']}
