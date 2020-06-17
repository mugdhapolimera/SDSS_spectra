# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 11:42:56 2020

@author: mugdhapolimera

Extract and plot MANGA data
https://www.sdss.org/dr16/manga/manga-tutorials/dap-tutorial/dap-python-tutorial/
www.sdss.org/dr16/manga/manga-data/data-model/
https://data.sdss.org/sas/dr16/manga/spectro/analysis/v2_4_3/2.2.1/HYB10-GAU-MILESHC/8655/6101/

"""
from astropy.io import fits
import os
import numpy as np
import matplotlib.pyplot as plt

os.chdir('F:\mugdhapolimera\Downloads')
f = fits.open('manga-8655-6101-LOGCUBE-HYB10-GAU-MILESHC.fits.gz')
hdu = fits.open('manga-8655-6101-MAPS-HYB10-GAU-MILESHC.fits.gz')
spec = f['emline'].data
wave = f['wave'].data

emline = {}
for k, v in hdu['EMLINE_GFLUX'].header.items():
    if k[0] == 'C':
        try:
            i = int(k[1:])-1
        except ValueError:
            continue
        emline[v] = i
print(emline) 
#plt.figure()
#plt.imshow(hdu['EMLINE_GFLUX'].data[emline['Ha-6564'],:,:],
#              cmap='inferno', origin='lower', interpolation='nearest')
#plt.colorbar(label=r'H$\alpha$ flux ($1\times10^{-17}$ erg s$^{-1}$ spaxel$^{-1}$ cm$^{-2}$)')
plt.figure()
mask_extension = hdu['EMLINE_GFLUX'].header['QUALDATA']
masked_halpha = np.ma.array(hdu['EMLINE_GFLUX'].data[emline['Ha-6564'],:,:],
                mask=hdu[mask_extension].data[emline['Ha-6564'],:,:]>0)
masked_oi = np.ma.array(hdu['EMLINE_GFLUX'].data[emline['OI-6365'],:,:],
                mask=hdu[mask_extension].data[emline['OI-6365'],:,:]>0)
plt.imshow(masked_halpha, origin='lower', cmap='inferno', interpolation='nearest')
plt.colorbar()
plt.show()

mask_ext = hdu['EMLINE_GVEL'].header['QUALDATA']
gas_vfield = np.ma.MaskedArray(hdu['EMLINE_GVEL'].data[emline['Ha-6564'],:,:], 
                    mask=hdu[mask_ext].data[emline['Ha-6564'],:,:] > 0)
#plt.figure()
#plt.imshow(gas_vfield, origin='lower', interpolation='nearest', vmin=-125, vmax=125, cmap='RdBu_r')
#plt.colorbar()
#plot the spectrum of the bin with the highest S/N
#snr = np.ma.MaskedArray(hdu['BIN_SNR'].data, mask=hdu['BINID'].data[0,:,:] < 0)
#j,i = np.unravel_index(snr.argmax(), snr.shape)
#
#flux = np.ma.MaskedArray(f['FLUX'].data[:,j,i],mask=f['MASK'].data[:,j,i] > 0)
#model = np.ma.MaskedArray(f['MODEL'].data[:,j,i],mask=f['MASK'].data[:,j,i]>0)
#stellarcontinuum = np.ma.MaskedArray(f['MODEL'].data[:,j,i] - f['EMLINE'].data[:,j,i] - f['EMLINE_BASE'].data[:,j,i], mask=f['MASK'].data[:,j,i] > 0)
#emlines = np.ma.MaskedArray(f['EMLINE'].data[:,j,i],mask=f['EMLINE_MASK'].data[:,j,i] > 0)
#resid = flux-model-0.5
#pyplot.clf()
#pyplot.step(wave, flux, where='mid', color='k', lw=0.5,label='flux')
#pyplot.plot(wave, model, color='r', lw=1,label='model')
#pyplot.plot(wave, stellarcontinuum, color='g', lw=1,label='stellar cont.')
#pyplot.plot(wave, emlines, color='b', lw=1,label='Emission lines')
#pyplot.step(wave, resid, where='mid', color='0.5', lw=0.5,label='residuals')
#pyplot.legend()

# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:29:16 2019

@author: mugdhapolimera

Reading and plotting SAMI IFU emission line measurements on a BPT plot

"""

import numpy as np
import pandas as pd
from astropy.table import Table as table
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib import colors
import os
from scipy.io import readsav

#define demarcation function: log_NII_HA vs. log_OIII_HB
def n2hacompmin(log_NII_HA): #composite minimum line from equation 1, Kewley 2006
    return 1.3 + (0.61 / (log_NII_HA - 0.05))
def n2hamain(log_NII_HA): #main line for NII/H-alpha from equation 5, Kewley 2006
    return 1.19 + (0.61 / (log_NII_HA - 0.47))
def s2hamain(log_SII_HA): #main line for SII/H-alpha from equation 2, Kewley 2006
    return 1.30 + (0.72 / (log_SII_HA - 0.32))
def s2halinseyf(log_SII_HA): #liner/seyfert divider for SII/H-alpha
    return 0.76 + 1.89*log_SII_HA
def o1hamain(log_OI_HA): #main line for OI/H-alpha from equation 3, Kewley 2006
    return 1.33 + (0.73 / (log_OI_HA + 0.59))
def o1halinseyf(log_OI_HA): #liner/seyfert divider for OI/H-alpha
    return 1.3 + 1.18*log_OI_HA
def o1hacrit(log_OI_HA): #boundary for OI/H-alpha
    return -0.59
def ratioerror(num,num_err,den, den_err):
    #err = (num/den) * np.sqrt((num_err/num)**2 + (den_err/den)**2)
    #err = ((num_err/num) + (den_err/den))/np.log(10)
    err_num2 = (num_err/(num*np.log(10)))**2
    err_den2 = (den_err/(den*np.log(10)))**2
    #err = np.sqrt((num_err/(num*np.log(10)))**2 + (den_err/(den*np.log(10)))**2)
    return np.sqrt(err_num2 + err_den2)
def logerr(err,number):
    log_err = 0.434*(err/number)
    return log_err

resdata = readsav('C:/Users/mugdhapolimera/github/SDSS_spectra/resolvecatalog.dat')
galname = 'rf0503'
df = pd.read_csv('C:\Users\mugdhapolimera\github\SDSS_spectra\RESOLVE_snr5_master.csv')
df.index = df.name
error = 0
halpha = hdu['EMLINE_GFLUX'].data[emline['Ha-6564'],:,:]
oi = hdu['EMLINE_GFLUX'].data[emline['OI-6302'],:,:]
hbeta = hdu['EMLINE_GFLUX'].data[emline['Hb-4862'],:,:]
oiii = hdu['EMLINE_GFLUX'].data[emline['OIII-5008'],:,:]
sii_1 = hdu['EMLINE_GFLUX'].data[emline['SII-6718'],:,:]
sii_2 = hdu['EMLINE_GFLUX'].data[emline['SII-6732'],:,:]
sii = sii_1 + sii_2
nii1 = hdu['EMLINE_GFLUX'].data[emline['NII-6549'],:,:]
nii2 = hdu['EMLINE_GFLUX'].data[emline['NII-6585'],:,:]
nii = (nii1 + nii2)*3.0/4

halpha_err = np.sqrt(1/hdu['EMLINE_GFLUX_IVAR'].data[emline['Ha-6564'],:,:])
oi_err = np.sqrt(1/hdu['EMLINE_GFLUX_IVAR'].data[emline['OI-6302'],:,:])
hbeta_err = np.sqrt(1/hdu['EMLINE_GFLUX_IVAR'].data[emline['Hb-4862'],:,:])
oiii_err = np.sqrt(1/hdu['EMLINE_GFLUX_IVAR'].data[emline['OIII-5008'],:,:])
sii_1_err = np.sqrt(1/hdu['EMLINE_GFLUX_IVAR'].data[emline['SII-6718'],:,:])
sii_2_err = np.sqrt(1/hdu['EMLINE_GFLUX_IVAR'].data[emline['SII-6732'],:,:])
sii_err = np.sqrt(sii_1_err**2 + sii_2_err**2)
nii1_err = np.sqrt(1/hdu['EMLINE_GFLUX_IVAR'].data[emline['NII-6549'],:,:])
nii2_err = np.sqrt(1/hdu['EMLINE_GFLUX_IVAR'].data[emline['NII-6585'],:,:])
nii_err = np.sqrt(nii1_err**2 + nii2_err**2)*3.0/4

#pixelscale = 0.5 #arcsec/pix
hdu2 = hdu['EMLINE_GFLUX'].header
pixelscale = hdu2['PC2_2'] #deg/pix
v = resdata['vhel'][np.where(resdata['name'] == galname)] #km/s
z = v/3e5
d = (v/70)*10**6 #Hubble's constant
pc = 2*np.pi*d*(pixelscale/360)/1000 #kpc/pix
#pc = pixelscale *3600 #in arcsec
MAX = np.max(halpha[~np.isnan(halpha)])
center = [[int(hdu2['CRPIX1'])],[int(hdu2['CRPIX2'])]]#np.where(halpha == MAX)
#center = [center[0][0], center[1][0]]
dx = np.arange(np.shape(halpha)[0])- center[0]
dy = np.arange(np.shape(halpha)[1])- center[1]
dxy = np.array(zip(dx,dy))
b_a = np.float(df.loc[galname].b_a)
r = np.zeros(halpha.shape)
for i in range(len(dx)):
    for j in range(len(dy)):
        r[i][j] = np.sqrt((pc*dx[i])**2 + (pc*-dy[j]/b_a)**2)
rndx = np.where(r<=2.05)

lam = wave 
mesh = np.meshgrid(np.arange(50), np.arange(50))
#To search for a particular wavelength
ndx = np.where(abs(lam-6685) == min(abs(lam-6685)))[0][0]

ha_cen = halpha[center]
hb_cen = hbeta[center]
nii_cen = nii[center]
sii_cen = sii[center]
oi_cen = oi[center]
oiii_cen = oiii[center]
ha_cen_err = halpha_err[center]
hb_cen_err = hbeta_err[center]
nii_cen_err = nii_err[center]
#sii_cen_err = sii_err[center]
oi_cen_err = oi_err[center]
oiii_cen_err = oiii_err[center]
center = [center[0][0],center[1][0]]

nans = [~np.isnan(halpha) & ~np.isnan(hbeta) & ~np.isnan(nii) & ~np.isnan(sii_1)
        & ~np.isnan(sii_2) & ~np.isnan(oi) & ~np.isnan(oiii)]
err = [(nii/nii_err > 5) & (sii_1/sii_1_err > 5) & (sii_2/sii_2_err > 5) & 
       (oi/oi_err > 5) & (oiii/oiii_err > 5) & (hbeta/hbeta_err > 5) & 
       (halpha/halpha_err > 5)]
sdss = (r <= 2.05)
good = (nans and err)# and sdss)
from copy import copy
cmap = copy(plt.cm.seismic)
cmap.set_bad('gray',0.8)

from astropy import wcs
from astropy.coordinates import SkyCoord as sky
cube_hdu = hdu['EMLINE_GFLUX'].header
w = wcs.WCS(cube_hdu)
w = w.dropaxis(2)
xy = np.meshgrid(np.arange(50), np.arange(50))
coords = w.all_pix2world(xy[0], xy[1],0)

image = np.ma.masked_where(~good[0],np.log10(oi/halpha))
goodimage = np.log10(masked_oi/masked_halpha)

#goodimage = np.ma.masked_where(((bluesnr < 5) & (redsnr < 5)),np.log10(oi/halpha))
#goodimage = np.ma.masked_where((r > 16),np.log10(oi/halpha))
#goodimage[:,0:15] = np.nan
#goodimage[:,34:] = np.nan
fig = plt.figure()
ax = plt.subplot(projection = w)
#ax.imshow(cube[ndx,:,:],norm = colors.Normalize(vmin = 0, vmax = 0.04439), 
#          cmap = 'inferno')
cax = ax.imshow(goodimage,
                norm = colors.Normalize(vmin = -1.4, vmax = -1.0), 
                cmap = 'seismic')
cax = ax.imshow(image,norm = colors.Normalize(vmin = -1.4, vmax = -1.0), 
          cmap = cmap)
ax.plot(28,28,'x', c = 'black')
fig.colorbar(cax,extend = 'min')
ax.set_xlabel('RA')
ax.set_ylabel('Dec')

#good = pix
halpha = halpha[good]
hbeta = hbeta[good]
nii = nii[good]
sii_1 = sii_1[good]
sii_2 = sii_2[good]
oi = oi[good]
oiii = oiii[good]
r = r[good]
halpha_err = halpha_err[good]
hbeta_err = hbeta_err[good]
nii_err = nii_err[good]
sii_1_err = sii_1_err[good]
sii_2_err = sii_2_err[good]
oi_err = oi_err[good]
oiii_err = oiii_err[good]
sii = sii[good]
sii_err = sii_err[good]

#good = (r <= 5.0)
#halpha = halpha[good]
#hbeta = hbeta[good]
#nii = nii[good]
#sii_1 = sii_1[good]
#sii_2 = sii_2[good]
#oi = oi[good]
#oiii = oiii[good]
#r = r[good]
#halpha_err = halpha_err[good]
#hbeta_err = hbeta_err[good]
#nii_err = nii_err[good]
#sii_1_err = sii_1_err[good]
#sii_2_err = sii_2_err[good]
#oi_err = oi_err[good]
#oiii_err = oiii_err[good]
#sii = sii[good]
#sii_err = sii_err[good]

#catid = [folder]*len(halpha)
#data = list(zip(catid, hbeta, hbeta_err,oiii, oiii_err, oi, oi_err,
#                halpha, halpha_err, nii, nii_err, sii_1, sii_1_err, sii_2, sii_2_err))
#names = ['CATID',
#         'h_beta_flux', 'h_beta_flux_err', 
#       'oiii_5007_flux', 'oiii_5007_flux_err',
#       'oi_6300_flux', 'oi_6300_flux_err', 
#       'h_alpha_flux','h_alpha_flux_err',
#       'nii_6584_flux', 'nii_6584_flux_err', 
#       'sii_6717_flux','sii_6717_flux_err',
#       'sii_6731_flux', 'sii_6731_flux_err']

#df = pd.DataFrame(data, columns = names)
#df.to_csv(folder+'.csv')

n2ha = np.log10(nii/halpha)
s2ha = np.log10(sii/halpha)
o1ha= np.log10(oi/halpha)
o3hb = np.log10(oiii/hbeta)
#n2ha_err = logerr(ratioerror(nii, nii_err, halpha, halpha_err),n2ha)
#s2ha_err = logerr(ratioerror(sii, sii_err, halpha, halpha_err),s2ha)
#o1ha_err = logerr(ratioerror(oi, oi_err, halpha, halpha_err),o1ha)
#o3hb_err = logerr(ratioerror(oiii, oiii_err, hbeta, hbeta_err),o3hb)
n2ha_err = ratioerror(nii, nii_err, halpha, halpha_err)
s2ha_err = ratioerror(sii, sii_err, halpha, halpha_err)
o1ha_err = ratioerror(oi, oi_err, halpha, halpha_err)
o3hb_err = ratioerror(oiii, oiii_err, hbeta, hbeta_err)

refn2ha = np.linspace(-3.0, 0.35)
refoiha = np.linspace(-2.5, -0.4)
refsiiha = np.linspace(-2, 0.32,100)
xlims = [-1.0,0.0]
ylims = [0,0.35]
cmap = plt.get_cmap('rainbow',int(np.max(r)/0.5))
fig,(ax1,ax2,ax3) = plt.subplots(1,3,sharey = True)#'NII Scatter Plot')
ax1.plot(refn2ha, n2hamain(refn2ha), 'k', 
                  label = 'ke01 Theoretical Maximum Starburst Line')
ax1.plot(refn2ha[refn2ha < 0], n2hacompmin(refn2ha[refn2ha < 0]),
                      'k-.', label = 'Ka03 Composite Line')
#ax1.plot(n2ha, o3hb, 'k.', alpha = 0.5, markersize = 5)#, label = 'Definite Star Forming')
ax1.set_xlabel(r"$\rm \log([NII]/H\alpha)$", fontsize = 22)
ax1.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)

#SII/OIII plot
#fig,ax2 = plt.subplots(322)#'SII Scatter Plot')
ax2.plot(refsiiha, s2hamain(refsiiha), 'k',  label = 'Ke01 Line')
ax2.plot(refsiiha[refsiiha > -0.32], s2halinseyf(refsiiha[refsiiha > -0.32]),
                  'k--', label = 'Liner/Seyfert Division')
#ax2.plot(s2ha, o3hb, 'k.', markersize = 5, \
#                    alpha = 0.5, label = 'SF')
#ax2.plot(np.log10(sii_cen/ha_cen), np.log10(oiii_cen/hb_cen), 'ro', alpha = 0.5, markersize = 5)#, label = 'Definite Star Forming')
ax2.set_xlabel(r"$\rm \log([SII]/H\alpha)$", fontsize = 22)
#ax2.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)

#OI/OIII plot
#fig ,ax3= plt.subplots(323)#'OI Scatter Plot')
ax3.plot(refoiha[refoiha < -0.7], o1hamain(refoiha[refoiha < -0.7]),
                  'k', label = 'Ke01 Theoretical Maximum Starburst Line')
ax3.plot(refoiha[refoiha < -0.7], o1hamain(refoiha[refoiha < -0.7]),
                  'k-.', label = 'Ka03 Composite Line')
#ax3.set_ylim(ylims)
ax3.plot(refoiha[refoiha > -1.13], o1halinseyf(refoiha[refoiha > -1.13]),'k--', 
         label = 'Ke06 Liner/Seyfert Division Line')
ax3.set_xlabel(r"$\rm \log([OI]/H\alpha)$", fontsize = 22)
#ax3.set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
#ax3.plot(o1ha, o3hb, 'k.', alpha = 0.5, 
#                    markersize = 5, label = 'SF')
#ax3.plot(np.log10(oi_cen/ha_cen), np.log10(oiii_cen/hb_cen), 'ro', alpha = 0.5, markersize = 5)#, label = 'Definite Star Forming')
error = 1
if error:
    ax1.errorbar(n2ha.flatten(), o3hb.flatten(), xerr = n2ha_err.flatten(),
                yerr = o3hb_err.flatten(),fmt = 'None', marker = 'None', 
                alpha = 0.5, mew = 0, label = 'SF-to-AGN',
                ecolor = 'k', zorder=0)
    ax2.errorbar(s2ha.flatten(), o3hb.flatten(), xerr = s2ha_err.flatten(),
                yerr = o3hb_err.flatten(), fmt = 'None', marker = 'None', c = r,
                alpha = 0.5, mew = 0, label = 'SF-to-AGN', ecolor = 'k',
                zorder=0)
    ax3.errorbar(o1ha.flatten(), o3hb.flatten(), xerr = o1ha_err.flatten(),
                yerr = o3hb_err.flatten(), fmt = 'None', marker = 'None', c = r,
                alpha = 0.5, mew = 0, label = 'SF-to-AGN', ecolor = 'k', 
                zorder=0)
cax = ax1.scatter(n2ha, o3hb,c = r, cmap = cmap)
cax = ax2.scatter(s2ha, o3hb,c = r, cmap = cmap)
cax = ax3.scatter(o1ha, o3hb,c = r, cmap = cmap)
fig.colorbar(cax,extend = 'min')
ax1.plot(np.log10(nii_cen/ha_cen), np.log10(oiii_cen/hb_cen), 'ko', \
         mfc = 'none', mew = 2, markersize = 15)#, label = 'Definite Star Forming')
ax2.plot(np.log10(sii_cen/ha_cen), np.log10(oiii_cen/hb_cen), 'ko', \
         mfc = 'none', mew = 2, markersize = 15)#, label = 'Definite Star Forming')
ax3.plot(np.log10(oi_cen/ha_cen), np.log10(oiii_cen/hb_cen), 'ko', \
         mfc = 'none', mew = 2, markersize = 15)#, label = 'Definite Star Forming')
ax1.scatter(np.log10(nii_cen/ha_cen), np.log10(oiii_cen/hb_cen),marker = 'o', 
            s = 175, c = [0.0], cmap = cmap)#, label = 'Definite Star Forming')
ax2.scatter(np.log10(sii_cen/ha_cen), np.log10(oiii_cen/hb_cen), marker = 'o', \
         s = 175, c = [0.0], cmap = cmap)#, label = 'Definite Star Forming')
ax3.scatter(np.log10(oi_cen/ha_cen), np.log10(oiii_cen/hb_cen), marker = 'o', \
         s = 175, c = [0.0], cmap = cmap)#, label = 'Definite Star Forming')

#ax3.errorbar(-0.9, 0.3, xerr = ratioerror(oi_cen,oi_cen/10,ha_cen,ha_cen/20),
#            yerr = ratioerror(oiii_cen,oiii_cen/20,hb_cen,hb_cen/10), 
#                marker = 'o', c = 'k')
#rs1105
df = pd.read_csv('C:\Users\mugdhapolimera\github\SDSS_spectra\RESOLVE_full_snr5.csv')
df.index = df.name
sdss_oi = df.oi_6300_flux.loc[galname]
sdss_oi_err = df.oi_6300_flux_err.loc[galname]
sdss_oiii =  df.oiii_5007_flux.loc[galname]
sdss_oiii_err = df.oiii_5007_flux_err.loc[galname]
sdss_ha = df.h_alpha_flux.loc[galname]
sdss_ha_err = df.h_alpha_flux_err.loc[galname]
sdss_hb = df.h_beta_flux.loc[galname]
sdss_hb_err = df.h_beta_flux_err.loc[galname]

#ax3.errorbar(-1.05, 0.3, xerr = ratioerror(oi_cen,oi_cen_err,ha_cen,ha_cen_err),
#                yerr = ratioerror(oiii_cen,oiii_cen_err,hb_cen,hb_cen_err), 
#                marker = 'o', c = 'k')
ax3.errorbar(np.log10(sdss_oi/sdss_ha), 
             np.log10(sdss_oiii/sdss_hb), xerr = ratioerror(sdss_oi,sdss_oi_err,sdss_ha,sdss_ha_err),
            yerr = ratioerror(sdss_oiii,sdss_oiii_err,sdss_hb,sdss_hb_err), 
                marker = 's', c = 'green')
#ax3.errorbar(-0.9, 0.3, xerr = ratioerror(sdss_oi,sdss_oi/23,sdss_ha,sdss_ha/43),
#            yerr = ratioerror(sdss_oiii,sdss_oiii/30,sdss_hb,sdss_hb/28), 
#                marker = 'o', c = 'k')

#ax3.errorbar(-1.05, 0.3, xerr = ratioerror(oi_cen,oi_cen_err,ha_cen,ha_cen_err),
#                yerr = ratioerror(oiii_cen,oiii_cen_err,hb_cen,hb_cen_err), 
#                marker = 'o', c = 'k')
#ax3.errorbar(-1.2, 0.3, xerr = ratioerror(sdss_oi,sdss_oi_err,sdss_ha,sdss_ha_err),
#            yerr = ratioerror(sdss_oiii,sdss_oiii_err,sdss_hb,sdss_hb_err), 
#                marker = 'o', c = 'k')
#ax3.errorbar(-0.9, 0.3, xerr = ratioerror(sdss_oi,sdss_oi/23,sdss_ha,sdss_ha/43),
#            yerr = ratioerror(sdss_oiii,sdss_oiii/30,sdss_hb,sdss_hb/28), 
#                marker = 'o', c = 'k')
ax3.set_xlim(-1.6,-0.8)
ax2.set_xlim(-0.5,-0.1)
ax1.set_xlim(xlims)
ax1.set_ylim(ylims)

