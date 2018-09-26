##deredden SDSS flux measurements for RESOLVE
##October 7, 2016

##############################
###############################


import pyfits 
from astropy.io import fits
from scipy import ndimage
from scipy.io.idl import readsav
from scipy.optimize import curve_fit
import numpy as np
from time import clock
import glob
import pandas as pd

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

from numpy import pi
from numpy.ma import median
from matplotlib import pyplot as plt
import os.path
import sys
import pdb
import pylab
pylab.ion()

#import and concatenate all the batches

#t1 = fits.open('coords/resolve_SDSS_one_SNR3.fits')
#t2 = fits.open('coords/resolve_SDSS_two_SNR3.fits')

#nrows1 = t1[1].data.shape[0]
#nrows2 = t2[1].data.shape[0]
#nrows = nrows1 + nrows2 

#hdu = fits.new_table(fits.ColDefs(t1[1].columns), nrows=nrows)
#for name in t1[1].columns.names:
#    hdu.data.field(name)[nrows1:]=t2[1].data.field(name)
#hdu.writeto('RESOLVE_SDSS_raw_SNR3.fits', clobber= True)


#open the file
hdulist = fits.open('RESOLVE_SDSS_raw_SNR3.fits')
#dwarf
#hdulist = fits.open('coords/resolve_SDSS_dwarf_SNR3.fits')
#extract the data
hdu_data = hdulist[1].data
#extract the header
hdu_headers = hdulist[1].header
#separate the data columns
galname = hdu_data.field(0)
oii_3726_flux = hdu_data.field(1)
oii_3726_flux_err = hdu_data.field(2)
oii_3729_flux = hdu_data.field(3)
oii_3729_flux_err = hdu_data.field(4)
neiii_3869_flux = hdu_data.field(5)
neiii_3869_flux_err = hdu_data.field(6)
h_delta_flux = hdu_data.field(7)
h_delta_flux_err = hdu_data.field(8)
h_gamma_flux = hdu_data.field(9)
h_gamma_flux_err = hdu_data.field(10)
oiii_4363_flux = hdu_data.field(11)
oiii_4363_flux_err = hdu_data.field(12)
h_beta_flux = hdu_data.field(13)
h_beta_flux_err = hdu_data.field(14)
oiii_4959_flux = hdu_data.field(15)
oiii_4959_flux_err = hdu_data.field(16)
oiii_5007_flux = hdu_data.field(17)
oiii_5007_flux_err = hdu_data.field(18)
hei_5876_flux = hdu_data.field(19)
hei_5876_flux_err = hdu_data.field(20)
oi_6300_flux = hdu_data.field(21)
oi_6300_flux_err = hdu_data.field(22)
nii_6548_flux = hdu_data.field(23)
nii_6548_flux_err = hdu_data.field(24)
h_alpha_flux = hdu_data.field(25)
h_alpha_flux_err = hdu_data.field(26)
nii_6584_flux = hdu_data.field(27)
nii_6584_flux_err = hdu_data.field(28)
sii_6717_flux = hdu_data.field(29)
sii_6717_flux_err = hdu_data.field(30)
sii_6731_flux = hdu_data.field(31)
sii_6731_flux_err = hdu_data.field(32)
ariii_7135_flux = hdu_data.field(33)
ariii_7135_flux_err = hdu_data.field(34)
heii_4685_flux = hdu_data.field(39)
heii_4685_flux_err = hdu_data.field(40)
oii_3726_flux_port = hdu_data.field(41)
oii_3726_flux_port_err = hdu_data.field(42)
########reddening correct###############
#intrinsic balmer decrement
balm_dec_exp = 2.86

#observed balmer decrement
balm_dec_obs = h_alpha_flux/h_beta_flux

#relation from balmer decrement to color excess
cee = 3.1*(np.log10(balm_dec_obs) - np.log10(balm_dec_exp))

#color excess for each galaxy
EBV_excess = cee*0.77

#speed of light
c_kms = 3.0*10**5.0

#extinction curve points
odonnell_dat = np.loadtxt('odonnell94mwextcurve.txt', dtype = "float")
#milky way lambda
mw_lam = odonnell_dat[:,0]
#factor to relate each mw lam and color excess
mw_A_overEBV = odonnell_dat[:,2]

###################REST FRAME###################
oii_3726_flux_ext = np.zeros(len(EBV_excess))
oii_3729_flux_ext = np.zeros(len(EBV_excess))
neiii_3869_flux_ext = np.zeros(len(EBV_excess))
h_delta_flux_ext = np.zeros(len(EBV_excess))
h_gamma_flux_ext = np.zeros(len(EBV_excess))
oiii_4363_flux_ext = np.zeros(len(EBV_excess))
h_beta_flux_ext = np.zeros(len(EBV_excess))
oiii_4959_flux_ext = np.zeros(len(EBV_excess))
oiii_5007_flux_ext = np.zeros(len(EBV_excess))
hei_5876_flux_ext = np.zeros(len(EBV_excess))
oi_6300_flux_ext = np.zeros(len(EBV_excess))
nii_6548_flux_ext = np.zeros(len(EBV_excess))
h_alpha_flux_ext  = np.zeros(len(EBV_excess))
nii_6584_flux_ext = np.zeros(len(EBV_excess))
sii_6717_flux_ext = np.zeros(len(EBV_excess))
sii_6731_flux_ext = np.zeros(len(EBV_excess))
ariii_7135_flux_ext = np.zeros(len(EBV_excess))
heii_4685_flux_ext = np.zeros(len(EBV_excess))
oii_3726_flux_port_ext = np.zeros(len(EBV_excess))


#find conversion array for each galaxy at each wavelength
mw_A_atmwlam  = np.zeros(len(EBV_excess))

for i in np.arange(len(galname)):
    mw_A_atmwlam = mw_A_overEBV*EBV_excess[i]
    #there is a mw_A_atmwlam for each galaxy, this index is the same for all
    mw_A_oii_3726_rest = mw_A_atmwlam[242]
    mw_A_oii_3729_rest = mw_A_atmwlam[243]
    mw_A_neiii_3869_rest = mw_A_atmwlam[290]
    mw_A_h_delta_rest = mw_A_atmwlam[367]
    mw_A_h_gamma_rest = mw_A_atmwlam[447]
    mw_A_oiii_4363_rest = mw_A_atmwlam[454]
    mw_A_h_beta_rest = mw_A_atmwlam[620]
    mw_A_oiii_4959_rest = mw_A_atmwlam[653]
    mw_A_oiii_5007_rest = mw_A_atmwlam[669]
    mw_A_hei_5876_rest = mw_A_atmwlam[959]
    mw_A_oi_6300_rest = mw_A_atmwlam[1100]
    mw_A_nii_6548_rest = mw_A_atmwlam[1183]
    mw_A_h_alpha_rest = mw_A_atmwlam[1187]
    mw_A_nii_6584_rest = mw_A_atmwlam[1195]
    mw_A_sii_6717_rest = mw_A_atmwlam[1239]
    mw_A_sii_6731_rest = mw_A_atmwlam[1244]
    mw_A_ariii_7135_rest = mw_A_atmwlam[1378]
    mw_A_heii_4685_rest = mw_A_atmwlam[562]
    mw_A_oii_3726_port_rest = mw_A_atmwlam[242]
    #find deextinction for obs lambda for each galaxy
    oii_3726_deext_rest = 10.0**(mw_A_oii_3726_rest/2.5)
    oii_3729_deext_rest = 10.0**(mw_A_oii_3729_rest/2.5)
    neiii_3869_deext_rest = 10.0**(mw_A_neiii_3869_rest/2.5)
    h_delta_deext_rest = 10.0**(mw_A_h_delta_rest/2.5)
    h_gamma_deext_rest = 10.0**(mw_A_h_gamma_rest/2.5)
    oiii_4363_deext_rest = 10.0**(mw_A_oiii_4363_rest/2.5)
    h_beta_deext_rest = 10.0**(mw_A_h_beta_rest/2.5)
    oiii_4959_deext_rest = 10.0**(mw_A_oiii_4959_rest/2.5)
    oiii_5007_deext_rest = 10.0**(mw_A_oiii_5007_rest/2.5)
    hei_5876_deext_rest = 10.0**(mw_A_hei_5876_rest/2.5)
    oi_6300_deext_rest = 10.0**(mw_A_oi_6300_rest/2.5)
    nii_6548_deext_rest = 10.0**(mw_A_nii_6548_rest/2.5)
    h_alpha_deext_rest = 10.0**(mw_A_h_alpha_rest/2.5)
    nii_6584_deext_rest = 10.0**(mw_A_nii_6584_rest/2.5)
    sii_6717_deext_rest = 10.0**(mw_A_sii_6717_rest/2.5)
    sii_6731_deext_rest = 10.0**(mw_A_sii_6731_rest/2.5)
    ariii_7135_deext_rest = 10.0**(mw_A_ariii_7135_rest/2.5)
    heii_4685_deext_rest= 10.0**(mw_A_heii_4685_rest/2.5)
    oii_3726_port_deext_rest = 10.0**(mw_A_oii_3726_port_rest/2.5)
    #output is dextincted avg flux for each line in each galaxy
    oii_3726_flux_ext[i] = oii_3726_flux[i]*oii_3726_deext_rest
    oii_3729_flux_ext[i] = oii_3729_flux[i]*oii_3729_deext_rest
    neiii_3869_flux_ext[i] = neiii_3869_flux[i]*neiii_3869_deext_rest
    h_delta_flux_ext[i] = h_delta_flux[i]*h_delta_deext_rest
    h_gamma_flux_ext[i] = h_gamma_flux[i]*h_gamma_deext_rest
    oiii_4363_flux_ext[i] = oiii_4363_flux[i]*oiii_4363_deext_rest
    h_beta_flux_ext[i] = h_beta_flux[i]*h_beta_deext_rest
    oiii_4959_flux_ext[i] = oiii_4959_flux[i]*oiii_4959_deext_rest
    oiii_5007_flux_ext[i] = oiii_5007_flux[i]*oiii_5007_deext_rest
    hei_5876_flux_ext[i] = hei_5876_flux[i]*hei_5876_deext_rest
    oi_6300_flux_ext[i] = oi_6300_flux[i]*oi_6300_deext_rest
    nii_6548_flux_ext[i] = nii_6548_flux[i]*nii_6548_deext_rest
    h_alpha_flux_ext[i]  = h_alpha_flux[i]*h_alpha_deext_rest
    nii_6584_flux_ext[i] = nii_6584_flux[i]*nii_6584_deext_rest
    sii_6717_flux_ext[i] = sii_6717_flux[i]*sii_6717_deext_rest
    sii_6731_flux_ext[i] = sii_6731_flux[i]*sii_6731_deext_rest
    ariii_7135_flux_ext[i] = ariii_7135_flux[i]*ariii_7135_deext_rest
    heii_4685_flux_ext[i] = heii_4685_flux[i]*heii_4685_deext_rest
    oii_3726_flux_port_ext[i] = oii_3726_flux_port[i]*oii_3726_port_deext_rest

#export to fits file
col1 = fits.Column(name = 'NAME', format = '20A', array=np.array(galname))
col2 = fits.Column(name = 'oii_3726_flux_ext', format = 'E', array=oii_3726_flux_ext)
col3 = fits.Column(name = 'oii_3726_flux_ext_err', format = 'E', array=oii_3726_flux_err)
col4 = fits.Column(name = 'oii_3729_flux_ext', format = 'E', array=oii_3729_flux_ext)
col5 = fits.Column(name = 'oii_3729_flux_ext_err', format = 'E', array=oii_3729_flux_err)
col6 = fits.Column(name = 'neiii_3869_flux_ext', format = 'E', array=neiii_3869_flux_ext)
col7 = fits.Column(name = 'neiii_3869_flux_ext_err', format = 'E', array=neiii_3869_flux_err)
col8 = fits.Column(name = 'h_delta_flux_ext', format = 'E', array=h_delta_flux_ext)
col9 = fits.Column(name = 'h_delta_flux_ext_err', format = 'E', array=h_delta_flux_err)
col10 = fits.Column(name = 'h_gamma_flux_ext', format = 'E', array=h_gamma_flux_ext)
col11 = fits.Column(name = 'h_gamma_flux_ext_err', format = 'E', array=h_gamma_flux_err)
col12 = fits.Column(name = 'oiii_4363_flux_ext', format = 'E', array=oiii_4363_flux_ext)
col13 = fits.Column(name = 'oiii_4363_flux_ext_err', format = 'E', array=oiii_4363_flux_err)
col14 = fits.Column(name = 'h_beta_flux_ext', format = 'E', array=h_beta_flux_ext)
col15 = fits.Column(name = 'h_beta_flux_ext_err', format = 'E', array=h_beta_flux_err)
col16 = fits.Column(name = 'oiii_4959_flux_ext', format = 'E', array=oiii_4959_flux_ext)
col17 = fits.Column(name = 'oiii_4959_flux_ext_err', format = 'E', array=oiii_4959_flux_err)
col18 = fits.Column(name = 'oiii_5007_flux_ext', format = 'E', array=oiii_5007_flux_ext)
col19 = fits.Column(name = 'oiii_5007_flux_ext_err', format = 'E', array=oiii_5007_flux_err)
col20 = fits.Column(name = 'hei_5876_flux_ext', format = 'E', array=hei_5876_flux_ext)
col21 = fits.Column(name = 'hei_5876_flux_ext_err', format = 'E', array=hei_5876_flux_err)
col22 = fits.Column(name = 'oi_6300_flux_ext', format = 'E', array=oi_6300_flux_ext)
col23 = fits.Column(name = 'oi_6300_flux_ext_err', format = 'E', array=oi_6300_flux_err)
col24 = fits.Column(name = 'nii_6548_flux_ext', format = 'E', array=nii_6548_flux_ext)
col25 = fits.Column(name = 'nii_6548_flux_ext_err', format = 'E', array=nii_6548_flux_err)
col26 = fits.Column(name = 'h_alpha_flux_ext', format = 'E', array=h_alpha_flux_ext)
col27 = fits.Column(name = 'h_alpha_flux_ext_err', format = 'E', array=h_alpha_flux_err)
col28 = fits.Column(name = 'nii_6584_flux_ext', format = 'E', array=nii_6584_flux_ext)
col29 = fits.Column(name = 'nii_6584_flux_ext_err', format = 'E', array=nii_6584_flux_err)
col30 = fits.Column(name = 'sii_6717_flux_ext', format = 'E', array=sii_6717_flux_ext)
col31 = fits.Column(name = 'sii_6717_flux_ext_err', format = 'E', array=sii_6717_flux_err)
col32 = fits.Column(name = 'sii_6731_flux_ext', format = 'E', array=sii_6731_flux_ext)
col33 = fits.Column(name = 'sii_6731_flux_ext_err', format = 'E', array=sii_6731_flux_err)
col34 = fits.Column(name = 'ariii_7135_flux_ext', format = 'E', array=ariii_7135_flux_ext)
col35 = fits.Column(name = 'ariii_7135_flux_ext_err', format = 'E', array=ariii_7135_flux_err)
col36 = fits.Column(name = 'Flux_HeII_4685_ext', format = 'E', array=heii_4685_flux_ext)
col37 = fits.Column(name = 'Flux_HeII_4685__ext_Err', format = 'E', array =heii_4685_flux_err)
col38 = fits.Column(name = 'Flux_OII_3726_ext', format = 'E', array=oii_3726_flux_port_ext)
col39 = fits.Column(name = 'Flux_OII_3726_ext_Err', format = 'E', array=oii_3726_flux_port_err)
col40 = fits.Column(name = 'balmer_decrement', format = 'E', array=balm_dec_obs)
col41 = fits.Column(name = 'balmer_decrement_ext', format = 'E', array=(h_alpha_flux_ext/h_beta_flux_ext))


#condense colummns, construct new table, write out results to binary fits
cols = fits.ColDefs([col1, col2,col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, col18, col19, col20, col21,col22,col23,col24,col25,col26,col27,col28,col29,col30,col31,col32,col33,col34,col35, col36, col37, col38, col39, col40, col41])
tbhdu = fits.new_table(cols)
#tbhdu.writeto('RESOLVE_SDSS_dext.fits', clobber=True)
#tbhdu.writeto('RESOLVE_SDSS_dext_dwarf_SNR3.fits', clobber=True)
#####################################################################



################################################

#plot check BPT
#define demarcation function, log_NII_HA vs. log_OIII_HB
def log_OIII_HB_NII(log_NII_HA):
    return 1.3 + (0.61 / (log_NII_HA - 0.05))

def comp_OIII_HB_NII(log_NII_HA):
    return 1.19 + 0.61 / (log_NII_HA - 0.47)

# create line ratios [NII]/Halpha and [OIII]/Hbeta
#check with Chris on this value, and Keweley 01 
#nii_sum = nii_6584_flux_ext
nii_sum = (nii_6548_flux_ext + nii_6584_flux_ext)*(3./4) 

x = np.log10(nii_sum/h_alpha_flux_ext)
y = np.log10(oiii_5007_flux_ext/h_beta_flux_ext)


#create starforming vs. AGN line
Measured_Predicted_OIII_HB = log_OIII_HB_NII(x)
Measured_Predicted_comp_OIII_HB = comp_OIII_HB_NII(x)

#generate list of points in x direction to plot demarcation line over entire range of plot
Predicted_NII_HA = np.linspace(-3.0, 0.35)
#evaluate function at those points
Predicted_log_OIII_HB_NII = log_OIII_HB_NII(Predicted_NII_HA)
Predicted_comp_log_OIII_HB_NII = comp_OIII_HB_NII(Predicted_NII_HA)

############################select out spheroids for these line ratio ranges

###single out E/S0s, and plot where there is available data
# create list of known names and turn into array
spheroids = np.loadtxt('blue_es0.txt', dtype = np.str)
spher_arr = np.array(spheroids)

#create list of  names of galaxies chosen by the HEII diagram
active_stars =np.loadtxt('HeII_gals.txt', dtype = np.str)
active_arr = np.array(active_stars)

#create list of names of galaxies chosen by BPT diagram
bpt_gals = np.loadtxt('BPT_gals.txt', dtype = np.str)
bpt_arr = np.array(bpt_gals)

x_sel = np.zeros(len(spher_arr))
y_sel = np.zeros(len(spher_arr))
activex = np.zeros(len(active_arr))
activey = np.zeros(len(active_arr))
bptx = np.zeros(len(bpt_arr))
bpty = np.zeros(len(bpt_arr))
for i in np.arange(len(spher_arr)):
    if spher_arr[i] in galname:
        sel_gal = spheroids[i]
        x_sel[i]= x[(np.where(galname == sel_gal))]
        y_sel[i]= y[(np.where(galname == sel_gal))]
for i in np.arange(len(active_arr)):
    if active_arr[i] in galname:
        active_sel = active_stars[i]
        activex[i] = x[(np.where(galname == active_sel))]
        activey[i] = y[(np.where(galname == active_sel))]
for i in np.arange(len(bpt_arr)):
    if bpt_arr[i] in galname:
        bpt_sel = bpt_gals[i]
        bptx[i] = x[(np.where(galname == bpt_sel))]
        bpty[i] = y[(np.where(galname == bpt_sel))]


### sample regions above and below line to plot the points above and below separately, ie in different colors 
#galaxies above line (AGN)
Above_Predicted_AGN = y > Measured_Predicted_OIII_HB
Above_Predicted_comp = y > Measured_Predicted_comp_OIII_HB


#select range to avoid ugly log function behavior at origin of plot 
sel = Predicted_NII_HA < 0

#remove spurious point at (0,0)
no_zero = np.where(x_sel != 0.)
no_zero_active = np.where(activex != 0.) 
no_zero_bpt = np.where(bptx != 0.)

#plot
fig = plt.figure(1)
plt.clf
ax = fig.add_subplot(111)
demarc, = ax.plot(Predicted_NII_HA[sel], Predicted_log_OIII_HB_NII[sel], color = "black",linestyle = '--')
composite, = ax.plot(Predicted_NII_HA, Predicted_comp_log_OIII_HB_NII, color = "black", linewidth = 1.5)
full, = ax.plot(x, y, 'k.', markersize =4)
es0s, = ax.plot(x_sel[no_zero],y_sel[no_zero],'bo', markersize = 8)
active, = ax.plot(activex[no_zero_active],activey[no_zero_active],"rs", markersize = 8)
bpt, = ax.plot(bptx[no_zero_bpt],bpty[no_zero_bpt],"c+", markersize = 10, mew=2)
ax.set_xlim(-3,2)
ax.set_ylim(-3,4)
ax.set_xlabel(r"$\rm \log [NII] \lambda6584 / H\alpha $", fontsize = 18)
ax.set_ylabel(r"$\rm \log [OIII]\lambda5007 / H\beta $", fontsize = 18)
#plt.legend([full, es0s,active], ['RESOLVE galaxies', 'ES0s', 'HeII AGN'], loc = 'lower left', numpoints =1)
#plt.savefig("SDSS_RESOLVE_BPT.eps")

#plot
#fig, ax = plt.subplots(1,2)
#demarc, = ax[0].plot(Predicted_NII_HA[sel], Predicted_log_OIII_HB_NII[sel], color = "black",linestyle = '--')
#composite, = ax[0].plot(Predicted_NII_HA, Predicted_comp_log_OIII_HB_NII, color = "black", linewidth = 1.5)
#full, = ax[0].plot(x, y, 'k.', markersize =4)
#es0s, = ax[0].plot(x_sel[no_zero],y_sel[no_zero],'bo', markersize = 4)
#active, = ax[0].plot(activex[no_zero_active],activey[no_zero_active],"ro", markersize = 8)
#ax[0].set_xlim(-3,2)
#ax[0].set_ylim(-3,4)
#ax[0].set_xlabel(r"$\rm \log([NII]/H\alpha)$", fontsize = 22)
#ax[0].set_ylabel(r"$\rm \log([OIII]/H\beta)$", fontsize = 22)
#ax[0].set_title("SDSS RESOLVE BPT", fontsize = 20)
#ax[0].set(aspect = 'equal')
plt.savefig("SDSS_RESOLVE_BPT.eps")


###################################################################
#plot HeII 

#define AGN demarcation lines
def log_HEII_HB(log_NII_HA):
    return  -1.22 + 1.0 / ((8.92*log_NII_HA) + 1.32)


heii_sel = np.where(heii_4685_flux_ext > heii_4685_flux_err*3.)
galname_heii = galname[heii_sel]

#define axis
x_heii = np.log10(nii_sum[heii_sel]/h_alpha_flux_ext[heii_sel])
y_heii = np.log10(heii_4685_flux_ext[heii_sel]/h_beta_flux_ext[heii_sel])


#create starforming vs. AGN line
Measured_Predicted_HEII_HB = log_HEII_HB(x_heii)


#generate list of points in x direction to plot demarcation line over entire range of plot
#evaluate function at those points
#generate list of points in x direction to plot demarcation line over entire range of plot
Predicted_NII_HA = np.linspace(-3.0, 0.35)
Predicted_log_HEII_HB = log_HEII_HB(Predicted_NII_HA)


#####select out spheroids for these line ratio ranges

###single out E/S0s, and plot where there is available data
# create list of known names and turn into array
spheroids = np.loadtxt('blue_es0.txt', dtype = np.str)
spher_arr = np.array(spheroids)

#create list of  names of galaxies chosen by the HEII diagram
#resolve
active_stars =np.loadtxt('HeII_gals.txt', dtype = np.str)
active_arr = np.array(active_stars)
x_sel_heii = np.zeros(len(spher_arr), dtype = object)
y_sel_heii = np.zeros(len(spher_arr), dtype = object)
activex_heii = np.zeros(len(active_arr), dtype = object)
activey_heii = np.zeros(len(active_arr), dtype = object)
bptx_heii = np.zeros(len(bpt_arr), dtype = object)
bpty_heii = np.zeros(len(bpt_arr), dtype = object)
for i in np.arange(len(spher_arr)):
    if spher_arr[i] in galname_heii:
        sel_gal = spheroids[i]
        x_sel_heii[i]= x_heii[(np.where(galname_heii == sel_gal))]
        y_sel_heii[i]= y_heii[(np.where(galname_heii == sel_gal))]
for i in np.arange(len(active_arr)):
    if active_arr[i] in galname_heii:
        active_sel = active_stars[i]
        activex_heii[i] = x_heii[(np.where(galname_heii == active_sel))]
        activey_heii[i] = y_heii[(np.where(galname_heii == active_sel))]
for i in np.arange(len(bpt_arr)):
    if bpt_arr[i] in galname_heii:
        bpt_sel = bpt_gals[i]
        bptx_heii[i] = x_heii[(np.where(galname_heii == bpt_sel))]
        bpty_heii[i] = y_heii[(np.where(galname_heii == bpt_sel))]

#select range to avoid ugly log function behavior at origin of plot 
sel = Predicted_NII_HA < -0.15

### sample regions above and below line to plot the points above and below separately, ie in different colors 
#galaxies above line (AGN)
Above_Predicted_AGN_heii = y_heii > Measured_Predicted_HEII_HB




#remove spurious point at (0,0)
no_zero_heii = np.where(x_sel_heii != 0) 
no_zero_active_heii = np.where(activex_heii != 0) 
no_zero_bpt_heii = np.where(bptx_heii != 0.)

#plot
fig = plt.figure(2)
plt.clf
ax_heii = fig.add_subplot(111)
demarc, = ax_heii.plot(Predicted_NII_HA[sel], Predicted_log_HEII_HB[sel], color = "black",linewidth = 1.5)
full, = ax_heii.plot(x_heii, y_heii, 'k.', markersize =4)
es0s, = ax_heii.plot(x_sel_heii[no_zero_heii],y_sel_heii[no_zero_heii],'bo', markersize = 8)
active, = ax_heii.plot(activex_heii[no_zero_active_heii],activey_heii[no_zero_active_heii],"rs", markersize = 8)
bpt, = ax_heii.plot(bptx_heii[no_zero_bpt_heii],bpty_heii[no_zero_bpt_heii],"c+", markersize = 10, mew=2)
ax_heii.set_xlim(-3.0,1.0)
ax_heii.set_ylim(-3.0,1.0)
ax_heii.set_xlabel(r"$\rm \log [NII]\lambda6584 / H\alpha$", fontsize = 18)
ax_heii.set_ylabel(r"$\rm \log [HeII]\lambda4686 / H\beta$", fontsize = 18)
plt.legend([full, es0s,active, bpt], ['RESOLVE galaxies', 'ES0s', 'HeII AGN', 'BPT AGN'], loc = 'lower left', numpoints =1, frameon = False, prop={'size':12})
plt.savefig("SDSS_RESOLVE_HeII.eps")

#demarc, = ax[1].plot(Predicted_NII_HA[sel], Predicted_log_HEII_HB[sel], color = "black",linewidth = 1.5)
#full, = ax[1].plot(x_heii, y_heii, 'k.', markersize =4)
#es0s, = ax[1].plot(x_sel_heii[no_zero_heii],y_sel_heii[no_zero_heii],'bo', markersize = 8)
#active, = ax[1].plot(activex_heii[no_zero_active_heii],activey_heii[no_zero_active_heii],"ro", markersize = 8)
#ax[1].set_xlim(-3.0,1.0)
#ax[1].set_ylim(-3.0,1.0)
#ax[1].set_xlabel(r"$\rm \log([NII]/H\alpha)$", fontsize = 22)
#ax[1].set_ylabel(r"$\rm \log([HeII]/H\beta)$", fontsize = 22)
#ax[1].set_title("SDSS RESOLVE HeII", fontsize = 20)
#ax[1].set(aspect = 'equal')
#plt.tight_layout()
#plt.savefig("SDSS_RESOLVE_HeII.eps")

