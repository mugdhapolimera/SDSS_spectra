# -*- coding: utf-8 -*-
"""
Created on Sun Dec 02 16:30:50 2018

SMC + Milky Way Blended Extinction Correction

This code blends SMC and MW based extinction correction, based on the 
galaxy mass. 

    - SMC extinction curve (Gordon03) for dwarf galaxies with M < 10^9
    - MW extinction curve (O’Donnell94) for giant galaxies with M > 10^10
    - A linear combination of SMC and MW extinction for galaxies with 
      10^9 < M < 10^10
                A = logM % 9 , B = 1-A
                Dereddened flux = B * SMC extinction + A * MW extinction

@author: mugdhapolimera
"""
#from astropy.table import Table
import numpy as np
import pandas as pd


def dext(smcfile, mwfile, outputfile):
    #smcfile = 'C:\Users\mugdhapolimera\github\SDSS_Spectra\NSA_RESOLVE_full_smcdext.csv'
    #smcfile = r'C:\Users\mugdhapolimera\github\xray\XMM_AGN_smcdext.pkl'
    #smcfile = r'C:\Users\mugdhapolimera\Desktop\UNC\Courses\Research\SAMI Data\71146\71146_smcdext.csv'
    #smc0 = Table.read(smcfile, format='fits')
    smc = pd.read_csv(smcfile)
    smc.index = smc.name
    #mwfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra\NSA_RESOLVE_full_smcdext.csv'
    #mwfile = r'C:\Users\mugdhapolimera\Desktop\UNC\Courses\Research\SAMI Data\71146\71146_mwdext.csv'
    #mw0 = Table.read(mwfile, format='fits')
    mw = pd.read_csv(mwfile)
    mw.index = mw.name
    #resolve = pd.read_csv('RESOLVE_liveOctober2018.csv')
    #resolve.index = resolve.name
    ext_corr = smc.copy()
    jhucolndx = [x for x in range(len(ext_corr.columns)) if 'flux' in list(ext_corr.columns)[x]]
    portcolndx = [x for x in range(len(ext_corr.columns)) if 'flux' in list(ext_corr.columns)[x]]
    for gal in smc.index.values:
        print(gal)    
        if smc.logmstar.loc[gal] < 9:
            ext_corr.loc[gal] = smc.loc[gal]
        elif smc.logmstar.loc[gal] > 10:
            ext_corr.loc[gal] = mw.loc[gal]
        else:
            A = smc.logmstar.loc[gal] % 9
            B = 1 - A
            #print 'Blending for Galaxy ', gal, A, B
            #print smc.loc[gal][1:73]
            #if ~np.isnan(np.array(ext_corr.loc[gal][96:172],dtype = float)).all():
            #    ext_corr.loc[gal][96:172] = A * smc.loc[gal][96:172] + B * mw.loc[gal][96:172]
            if ~np.isnan(np.array(ext_corr.loc[gal][jhucolndx],dtype = float)).all(): #[150:225]
                ext_corr.loc[gal][jhucolndx] = A * smc.loc[gal][jhucolndx] + B * mw.loc[gal][jhucolndx]#[150:225]
            if ~np.isnan(np.array(ext_corr.loc[gal][portcolndx],dtype = float)).all(): #[150:225]
                ext_corr.loc[gal][portcolndx] = A * smc.loc[gal][portcolndx] + B * mw.loc[gal][portcolndx]#[150:225]
    print(ext_corr)
    #ext_corr.to_pickle('XMM_AGN_blend.pkl')
    #ext_corr.to_csv('XMM_AGN_blend.csv')
    #ext_corr.to_csv('C:\Users\mugdhapolimera\Desktop\UNC\Courses\Research\SAMI Data\71146\71147_blend.csv')
    ext_corr.to_csv(outputfile)#'NSA_RESOLVE_full_blend_dext.csv')
    #ext_corr.index = np.array(ext_corr['name'])
    #dfres = pd.merge(resolve,ext_corr,how="left",left_index=True,right_index=True)
    return ext_corr