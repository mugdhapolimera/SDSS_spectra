# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 03:09:33 2022

@author: mugdhapolimera
"""

import pandas as pd
import numpy as np
from scipy.io import readsav

import matplotlib.pyplot as plt

eco = 0
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
salimsfr = salim.logsfr_sed

ressel= pd.read_csv("RESOLVE_snr5_dext_catsfrs.csv")
rescatfull = pd.read_csv(survey+"_inobssample.csv")
rescatfull.index = rescatfull.name
rescat = rescatfull#.loc[salim.name]
res = readsav('resolvecatalog_021622.dat')
resphot = readsav('resolvecatalogphot_021622.dat')
common, catndx, resndx = np.intersect1d(rescat.name, res.name, return_indices = True)
rescat['sfr_nuv'].iloc[catndx] = res['sfr_nuv'][resndx]
rescat['sfr_nuv_err'] = np.zeros(len(rescat))
rescat['sfr_nuv_err'].iloc[catndx] = np.array(res['sfr_nuv_err'][resndx]).byteswap().newbyteorder()
rescat['sfr_nuv_wise'] = np.zeros(len(rescat))
rescat['sfr_nuv_wise'].iloc[catndx] = np.array(res['sfr_nuv_wise'][resndx]).byteswap().newbyteorder()
rescat['sfr_nuv_wise_err'] = np.zeros(len(rescat))
rescat['sfr_nuv_wise_err'].iloc[catndx] = np.array(res['sfr_nuv_wise_err'][resndx]).byteswap().newbyteorder()
rescat['mw1w'] = np.zeros(len(rescat))
rescat['mw1w'].iloc[catndx] = np.array(resphot['mw1w'][resndx]).byteswap().newbyteorder()
rescat['emw1w'] = np.zeros(len(rescat))
rescat['emw1w'].iloc[catndx] = np.array(resphot['emw1w'][resndx]).byteswap().newbyteorder()

ressfr = np.log10(rescat.sfr_nuv.loc[salim.name])
#ressfr = rescat.sfr_nuv.loc[salim.name]
residual = (salimsfr - ressfr)
#err_residual = salim.err_logsfr_sed#*(10**salimsfr)/0.434/ressfr
ressfr_err = 0.434*rescat['sfr_nuv_err'].loc[salim.name]/rescat.sfr_nuv.loc[salim.name]    
err_residual = (salim.err_logsfr_sed**2 + ressfr_err**2)**0.5

plt.figure()
plt.errorbar(ressfr, residual, yerr = err_residual, fmt = 'o')
plt.xlabel('log10('+survey+' SFR_NUV)')
plt.ylabel('log10(Salim SFR_UV_OPTICAL_SED) - log10('+survey+' SFR_NUV)')
plt.axhline(y=0, c = 'k')
plt.axhline(y=np.nanmedian(residual), c = 'r', ls = 'dashed')
#trendout = np.array(rescat.name[(residual+err_residual < 0) & (residual < -0.5) & 
#                                (ressfr < -0.5)]) 
trendout = np.array(rescat.name[(residual+err_residual < 0) & (residual < -0.7) | 
                                (residual < 0) & (ressfr < -1.25)]) 
outliers = [x for x in trendout if x not in list(ressel.name)]
plt.errorbar(ressfr.loc[outliers], residual.loc[outliers], 
             yerr = err_residual.loc[outliers], fmt = 'o', color = 'black')

nonoutliers= np.setdiff1d(rescat.name, outliers)
plt.figure()
plt.errorbar(ressfr.loc[nonoutliers], residual.loc[nonoutliers], 
             yerr = err_residual.loc[nonoutliers], fmt = 'o')
plt.xlabel('log10('+survey+' SFR_NUV)')
plt.ylabel('log10(Salim SFR_UV_OPTICAL_SED) - log10('+survey+' SFR_NUV)')
plt.axhline(y=0, c = 'k')
plt.axhline(y=np.nanmedian(residual.loc[nonoutliers]), c = 'r', ls = 'dashed')
#Covariance in y direction
#xarrow = np.std(salimsfr - ressfr.loc[salim.name])
#yarrow = np.std(salimsfr - ressfr.loc[salim.name])
#diag = np.sqrt(xarrow**2 + yarrow**2)
#plt.arrow(-0.75,0,xarrow, -yarrow, zorder = 5)
##Covariance in x direction
xarrow = yarrow = np.std(ressfr_err.loc[salim.name])
diag = np.sqrt(xarrow**2 + yarrow**2)
plt.arrow(1,-2.5,xarrow, -yarrow, zorder = 5, length_includes_head = True,
          linewidth = 3)


residual = salimsfr
err_residual = salim.err_logsfr_sed
plt.figure()
plt.errorbar(ressfr, residual, yerr = err_residual, fmt = 'o')
plt.xlabel('log10('+survey+' SFR_NUV)')
plt.ylabel('log10(Salim SFR_UV_OPTICAL_SED)')
plt.errorbar(ressfr.loc[outliers], residual.loc[outliers], 
             yerr = err_residual.loc[outliers], fmt = 'o', color = 'black')
xaxis = np.arange(np.nanmin(ressfr) - 0.2, np.nanmax(ressfr)+0.2, 0.2)
plt.plot(xaxis,xaxis)


###############################################################################
#SSFR trends
resssfr = - rescatfull.logmstar + np.log10(rescatfull.sfr_nuv)
resssfr = resssfr.loc[salim.name]
salimssfr = - salim.col10 + salim.logsfr_sed
residual = (salimssfr - resssfr)
err_residual = (salim.err_logsfr_sed**2 + salim.col11**2)**0.5
plt.figure()
plt.errorbar(resssfr, residual, yerr = err_residual, fmt = 'o')
plt.xlabel('log10('+survey+' SSFR from NUV)')
plt.ylabel('log10(Salim SSFR from UV/OPTICAL SED) - log10('+survey+' SSFR from NUV)')
plt.axhline(y=0, c = 'k')
plt.axhline(y=np.nanmedian(residual), c = 'r', ls = 'dashed')
#plt.ylim(-2.5, 1.5)
plt.errorbar(resssfr.loc[outliers], residual.loc[outliers], 
             yerr = err_residual.loc[outliers], fmt = 'o', color = 'black')
xaxis = np.arange(np.nanmin(resssfr) - 0.2, np.nanmax(resssfr)+0.2, 0.2)
cov = -3 + xaxis*-(np.std(salimssfr - resssfr.loc[salim.name]))**2
plt.plot(xaxis, cov, zorder = 5, lw = 3)


resssfr = - rescatfull.logmstar + np.log10(rescatfull.sfr_nuv)
salimssfr = - salim.col10 + salim.logsfr_sed
plt.figure()
plt.subplot(1,2,1)
plt.scatter(rescatfull.logmstar, np.log10(rescatfull.sfr_nuv))
plt.scatter(rescatfull.logmstar.loc[outliers], 
            np.log10(rescatfull.sfr_nuv.loc[outliers]), c = 'none',
            edgecolor = 'k', linewidth = 3, marker = 's', s = 100)
plt.ylim(-3,1.5)
plt.xlabel('log(M*)')
plt.ylabel('log(RESOLVE SFR_NUV)')
plt.subplot(1,2,2)
plt.scatter(salim.col10, salim.logsfr_sed)
plt.scatter(salim.col10[outliers], salim.logsfr_sed.loc[outliers], c = 'none',
            edgecolor = 'k', linewidth = 3, marker = 's', s = 100)
plt.xlabel('Salim log(M*)')
plt.ylabel('Salim log(SFR SED)')
plt.ylim(-3,1.5)

plt.figure()
plt.subplot(1,2,1)
plt.scatter(rescatfull.logmstar, resssfr)
plt.scatter(rescatfull.logmstar.loc[outliers], 
            resssfr.loc[outliers], c = 'none',
            edgecolor = 'k', linewidth = 3, marker = 's', s = 100)
plt.ylim(-13.5,-8)
plt.xlabel('RESOLVE log(M*)')
plt.ylabel('log(RESOLVE SFR_NUV/M*)')
plt.subplot(1,2,2)
plt.scatter(salim.col10, salimssfr)
plt.scatter(salim.col10[outliers], salimssfr.loc[outliers], c = 'none',
            edgecolor = 'k', linewidth = 3, marker = 's', s = 100)
plt.xlabel('Salim log(M*)')
plt.ylabel('Salim log(SFR SED/M*)')
plt.ylim(-13.5,-8)

#plt.figure()
#plt.scatter(rescatfull.logmgas, np.log10(rescatfull.sfr_nuv))
#plt.scatter(rescatfull.logmgas.loc[outliers], 
#            np.log10(rescatfull.sfr_nuv.loc[outliers]))
#plt.xlabel('log(M_gas)')
#plt.ylabel('log(RESOLVE SFR_NUV)')

reslive = pd.read_csv("RESOLVE_live02Feb2022.csv")
reslive.index = reslive.name
reslive = reslive.loc[rescatfull.name]
resolve_sf = (reslive['grpcz']/reslive['cz'])**2.0
reslive['mhidet'] *= resolve_sf
reslive['emhidet'] *= resolve_sf
reslive['mhilim'] *= resolve_sf
reslive['loggastostar'] = np.log10(reslive.mhidet) - reslive.logmstar
reslive['loggastostar_lim'] = np.log10(reslive.mhilim) - reslive.logmstar
detections = reslive.mhidet>0
uplim = (reslive['mhilim']>0)

plt.figure()
plt.subplot(1,2,1)
plt.scatter(np.log10(reslive.mhidet.loc[detections]), 
            np.log10(reslive.sfr_nuv.loc[detections]), label = 'HI detections')
plt.scatter(np.log10(reslive.mhilim.loc[uplim]), 
            np.log10(reslive.sfr_nuv.loc[uplim]), marker = 'v', color = 'red',
            label = 'HI upper limits')
plt.legend(loc = 'upper left')
plt.scatter(np.log10(reslive.mhidet.loc[outliers]), 
            np.log10(reslive.sfr_nuv.loc[outliers]), marker = 's', c = 'none',
            s= 100, edgecolor = 'k', linewidth = 3)
plt.scatter(np.log10(reslive.mhilim.loc[outliers]), 
            np.log10(reslive.sfr_nuv.loc[outliers]), marker = 's', c = 'none',
            s= 100, edgecolor = 'k', linewidth = 3)
plt.xlabel('log(M_HI)')
plt.ylabel('log(RESOLVE SFR_NUV)')
plt.ylim(-2.5,1.5)

plt.subplot(1,2,2)
plt.scatter(reslive.loggastostar.loc[detections], 
            np.log10(reslive.sfr_nuv.loc[detections]))
plt.scatter(reslive.loggastostar_lim.loc[uplim], 
            np.log10(reslive.sfr_nuv.loc[uplim]), marker = 'v', color = 'red')
plt.scatter(reslive.loggastostar.loc[outliers], 
            np.log10(reslive.sfr_nuv.loc[outliers]), marker = 's', c = 'none',
            s= 100, edgecolor = 'k', linewidth = 3)
plt.scatter(reslive.loggastostar_lim.loc[outliers], 
            np.log10(reslive.sfr_nuv.loc[outliers]), marker = 's', c = 'none',
            s= 100, edgecolor = 'k', linewidth = 3)
plt.xlabel('log(M_HI/M*)')
plt.ylabel('log(RESOLVE SFR_NUV)')
plt.ylim(-2.5,1.5)

detections = reslive.name.loc[salim.name][reslive.loc[salim.name]['mhidet']>0]
uplim = reslive.name.loc[salim.name][reslive.loc[salim.name]['mhilim']>0]
plt.figure()
plt.subplot(1,2,1)
plt.scatter(np.log10(reslive.mhidet.loc[detections]), 
            salim.logsfr_sed.loc[detections], label = 'HI detections')
plt.scatter(np.log10(reslive.mhilim.loc[uplim]), 
            salim.logsfr_sed.loc[uplim], marker = 'v', color = 'red',
            label = 'HI upper limits')
plt.legend(loc = 'upper left')
plt.scatter(np.log10(reslive.mhidet.loc[outliers]), 
            salim.logsfr_sed.loc[outliers], marker = 's', c = 'none',
            s= 100, edgecolor = 'k', linewidth = 3)
plt.scatter(np.log10(reslive.mhilim.loc[outliers]), 
            salim.logsfr_sed.loc[outliers], marker = 's', c = 'none',
            s= 100, edgecolor = 'k', linewidth = 3)
plt.xlabel('log(M_HI)')
plt.ylabel('log(Salim SFR_SED)')
plt.ylim(-3.5,1.0)

plt.subplot(1,2,2)
plt.scatter(reslive.loggastostar.loc[detections], 
            salim.logsfr_sed.loc[detections])
plt.scatter(reslive.loggastostar_lim.loc[uplim], 
            salim.logsfr_sed.loc[uplim], marker = 'v', color = 'red')
plt.scatter(reslive.loggastostar.loc[outliers], 
            salim.logsfr_sed.loc[outliers], marker = 's', c = 'none',
            s= 100, edgecolor = 'k', linewidth = 3)
plt.scatter(reslive.loggastostar_lim.loc[outliers], 
            salim.logsfr_sed.loc[outliers], marker = 's', c = 'none',
            s= 100, edgecolor = 'k', linewidth = 3)
plt.xlabel('log(M_HI/M*)')
plt.ylabel('log(Salim SFR_SED)')
plt.ylim(-3.5,1.0)
###############################################################################
#Compare extinction
#resdatmags = pd.read_csv("magsandsfr_res_rereduced_withfuv.csv")
#resdatmags.index = resdatmags.name
#
#resdatext = -1*resdatmags.inextfuv.loc[salim.name]
#salimext = salim.col14
#plt.figure()
#plt.scatter(resdatext, (salimext - resdatext))
#plt.axhline()
#plt.xlabel('RESOLVE Re-reduced data FUV Extinction = smoothrestmag_fuv - deextrestmag_fuv')
#plt.ylabel('SALIM FUV ext - RESOLVE Re-reduced data FUV Ext')
#
#rescatmags = pd.read_csv("magsandsfr_res_catalogdata_withfuv.csv")
#rescatmags.index = rescatmags.name
#rescatext = -1*rescatmags.inextfuv.loc[salim.name]
#salimext = salim.col14
#plt.figure()
#plt.scatter(rescatext, (salimext - rescatext))
#plt.axhline()
#plt.xlabel('RESOLVE catalog data FUV Extinction = smoothrestmag_fuv - deextrestmag_fuv')
#plt.ylabel('SALIM FUV ext - RESOLVE catalog data FUV Ext')

#salim = pd.read_csv("RESOLVE-inobssample_GSWLC.csv")
#salim.index = salim.name
#salim = salim[(salim.log_sfr_allwise != -99) & (salim.flag_sed == 0)]
##salimsfr = np.log10(0.17*10**salim.log_sfr_allwise + 10**salim.logsfr_sed)
#salimsfr = salim.log_sfr_allwise
#ressfr = np.log10(rescat.sfr_nuv_wise.loc[salim.name] - rescat.sfr_nuv.loc[salim.name])
#ressfr_err = np.sqrt(rescat['sfr_nuv_wise_err'].loc[salim.name]**2 + rescat['sfr_nuv_err'].loc[salim.name]**2)
#ressfr_err = 0.434*ressfr_err/ressfr    
#err_residual = ressfr_err
#residual = (salimsfr - ressfr)
#plt.figure()
#plt.errorbar(ressfr, residual, yerr = err_residual, fmt = 'o')
#plt.xlabel('log10(RESOLVE SFR_NUV_WISE)')
#plt.ylabel('log10(SDSS SFR_WISE) - log10(RESOLVE SFR_NUV_WISE)')
#plt.axhline(y=0, c = 'k')
#plt.axhline(y=np.nanmedian(residual), c = 'r', ls = 'dashed')


resssfr = np.log10(rescat.sfr_nuv_wise * 1e9) - rescat.logmstar
 