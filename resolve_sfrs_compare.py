# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 08:20:35 2022

@author: mugdhapolimera
"""

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

#ressel= pd.read_csv("RESOLVE_snr5_dext_catsfrs.csv")
rescat = pd.read_csv(survey+"_live02Feb2022.csv")
rescat.index = rescat.name
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
rescat['mw2w'] = np.zeros(len(rescat))
rescat['mw2w'].iloc[catndx] = np.array(resphot['mw2w'][resndx]).byteswap().newbyteorder()
rescat['emw2w'] = np.zeros(len(rescat))
rescat['emw2w'].iloc[catndx] = np.array(resphot['emw2w'][resndx]).byteswap().newbyteorder()
rescat['mw3w'] = np.zeros(len(rescat))
rescat['mw3w'].iloc[catndx] = np.array(resphot['mw3w'][resndx]).byteswap().newbyteorder()
rescat['emw3w'] = np.zeros(len(rescat))
rescat['emw3w'].iloc[catndx] = np.array(resphot['emw3w'][resndx]).byteswap().newbyteorder()
rescat['mw4w'] = np.zeros(len(rescat))
rescat['mw4w'].iloc[catndx] = np.array(resphot['mw4w'][resndx]).byteswap().newbyteorder()
rescat['emw4w'] = np.zeros(len(rescat))
rescat['emw4w'].iloc[catndx] = np.array(resphot['emw4w'][resndx]).byteswap().newbyteorder()

#rescat = rescat[(rescat.sfr_nuv_wise_err < 100)]
sel = ~np.isnan(rescat.sfr_nuv_wise_err)
sfrnuv = np.log10(rescat.sfr_nuv.loc[sel])
sfrnuvwise = np.log10(rescat.sfr_nuv_wise.loc[sel])
residual = (sfrnuvwise-sfrnuv)
sfrnuv_err= (0.434*rescat['sfr_nuv_err']/rescat.sfr_nuv).loc[sel]    
sfrnuvwise_err= (0.434*rescat['sfr_nuv_wise_err']/rescat.sfr_nuv_wise).loc[sel]    
err_residual = (sfrnuv_err**2 + sfrnuvwise_err**2)**0.5

plt.figure()
#plt.errorbar(sfrnuv, residual, yerr = err_residual, fmt = 'o')
plt.scatter(sfrnuv, residual)
plt.xlabel('log10('+survey+' SFR_NUV)')
plt.ylabel('log10('+survey+' SFR_NUV_WISE) - log10('+survey+' SFR_NUV)')
plt.axhline(y=0, c = 'k')
plt.axhline(y=np.nanmedian(residual), c = 'r', ls = 'dashed')
#outliers = np.array(rescat.loc[sel].name[(residual> 0) & (residual-err_residual > 0)]) 
outliers = np.array(rescat.loc[sel].name[(residual> 0.3) & (residual-err_residual > 0)]) 
plt.errorbar(sfrnuv.loc[outliers], residual.loc[outliers], 
             yerr = err_residual.loc[outliers], fmt = 'o', color = 'red')


import scipy
from scipy.stats import kde
import matplotlib
import matplotlib.cm as cm
import matplotlib.colors as mpcolors

def density_estimation(m1, m2):
    X, Y = np.mgrid[xmin:xmax:25j, ymin:ymax:25j]                                                     
    positions = np.vstack([X.ravel(), Y.ravel()])                                                       
    values = np.vstack([m1, m2])                                                                        
    kernel = scipy.stats.gaussian_kde(values)                                                                 
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z 
plt.figure()
ax1 = plt.subplot()
ymin, ymax = (0, 3)
xmin,xmax = (7.0,11.5)

sel = (rescat.umag > 0) & (rescat.rmag > 0)
u_r = rescat.umag[sel] - rescat.rmag[sel]
mass = rescat.logmstar[sel]
X,Y,Z = density_estimation(mass,u_r)
Z = Z/Z.max()
lvls = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
def truncate_colormap(cmap, minval=0, maxval=0.75, n=150):
  	new_cmap = mpcolors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
  	return new_cmap
sf_colors_map = truncate_colormap(cm.gray,n=11)
nbins = 20
ax1.contour(X, Y, Z, levels = lvls, cmap=sf_colors_map, zorder = 4)
ax1.scatter(mass.loc[outliers], 
            u_r.loc[outliers], c = 'r')
ax1.set_xlabel('log(M*)')
ax1.set_ylabel('umag - rmag (no extinction correction)')
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
