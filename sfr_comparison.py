# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 16:53:39 2021

@author: mugdhapolimera
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.io.idl import readsav
eco = 0
sfr = 1
sel = 0
if eco:
    survey = 'ECO'
else:
    survey = 'RESOLVE'

salim = pd.read_csv(survey+"-inobssample_GSWLC.csv")
if not(eco) and sel:
    salim = pd.read_csv("RESOLVE-SEL_GSWLC.csv")
if eco and sel:
    salim = pd.read_csv("ECO-SEL_GSWLC.csv")

salim.index = salim.name
salim = salim[(salim.flag_sed == 0) & (salim.logsfr_sed != -99)]
#rescat = pd.read_csv("RESOLVE_snr5_dext_catsfrs.csv")
rescat = pd.read_csv(survey+"_inobssample.csv")
rescat.index = rescat.name
res = readsav('resolvecatalog_021222.dat')
common, catndx, resndx = np.intersect1d(rescat.name, res.name, return_indices = True)
rescat['sfr_nuv'].iloc[catndx] = res['sfr_nuv'][resndx]
rescat['sfr_nuv_err'] = np.zeros(len(rescat))
rescat['sfr_nuv_err'].iloc[catndx] = np.array(res['sfr_nuv_err'][resndx]).byteswap().newbyteorder()
##rescat = rescat[(10**rescat.logmgas + 10**rescat.logmstar) > 10**9.2]

if sfr:
    salimsfr = salim.logsfr_sed
    ressfr = np.log10(rescat.loc[salim.name].sfr_nuv)
    residual = (salimsfr - ressfr)
    if 'sfr_nuv_err' in rescat.keys():
        ressfr_err = 0.434*rescat['sfr_nuv_err'].loc[salim.name]/rescat.sfr_nuv.loc[salim.name]    
        err_residual = (salim.err_logsfr_sed**2 + ressfr_err**2)**0.5
                    #*(10**salimsfr)/0.434/ressfr
    else:
        err_residual = salim.err_logsfr_sed

    plt.figure()
    plt.errorbar(ressfr, residual, yerr = err_residual, fmt = 'o')
    plt.xlabel('log10('+survey+' SFR_NUV)')
    plt.ylabel('log10(SDSS SFR_UV_OPTICAL_SED) - log10('+survey+' SFR_NUV)')
    plt.axhline(y=0, c = 'k')
    plt.axhline(y=np.nanmedian(residual), c = 'r', ls = 'dashed')
else:
    salimsfr = salim.logsfr_sed - salim.col10
    ressfr = np.log10(rescat.sfr_nuv.loc[salim.name]) - rescat.logmstar.loc[salim.name]
    #ressfr = rescat.sfr_nuv.loc[salim.name]
    residual = (salimsfr - ressfr)
    err_residual = (salim.err_logsfr_sed**2 + salim.col11**2)**0.5 

    plt.figure()
    plt.errorbar(ressfr, residual, yerr = err_residual, fmt = 'o')
    plt.xlabel('log10('+survey+' sSFR_NUV)')
    plt.ylabel('log10(SDSS sSFR_UV_OPTICAL_SED) - log10('+survey+' sSFR_NUV)')
    plt.axhline(y=0, c = 'k')
    plt.axhline(y=np.nanmedian(residual), c = 'r', ls = 'dashed')


#d = rescat.cz/70.0*1e6*3.0857e18 #in cm
#L_ha = rescat.h_alpha_flux*1e-17 * 4 * np.pi * d**2
#sfr_ha = np.log10(7.9*1e-42*L_ha) - 0.24 
##sfr(Ha) from kennicutt98; 0.24dex offset to convert from salpeter to chabrier IMF (salim+16)
#ressfr = np.log10(rescat.sfr_nuv_wise)
#residual = (sfr_ha - ressfr)
##err_residual = salim.err_logsfr_sed
#plt.figure()
#plt.plot(ressfr, residual, 'o')
#plt.xlabel('log10(RESOLVE SFR_NUV_WISE)')
#plt.ylabel('log10(RESOLVE SFR(Ha)) - log10(RESOLVE SFR_NUV_WISE)')
#plt.axhline(y=0, c = 'k')
#plt.axhline(y=np.nanmedian(residual), c = 'r', ls = 'dashed')
#
#plt.figure()
#plt.hist(sfr_ha)
#plt.axvline(x = np.log10(0.0221), c = 'k', ls = 'dashed')
