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
#pylab.ion()

os.chdir('C:\Users\mugdhapolimera\github\SDSS_Spectra')
#open the file
hdulist = fits.open('ECO_full.fits')

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
#oiii_flux = hdu_data.field(39)
#oiii_flux_err = hdu_data.field(40)

heii_3203_flux_port = hdu_data.field(39)
heii_3203_flux_port_err = hdu_data.field(40)
oii_3726_flux_port = hdu_data.field(41)
oii_3726_flux_port_err = hdu_data.field(42)
oii_3729_flux_port = hdu_data.field(43)
oii_3729_flux_port_err = hdu_data.field(44)
neiii_3869_flux_port = hdu_data.field(45)
neiii_3869_flux_port_err = hdu_data.field(46)
h_delta_flux_port = hdu_data.field(47)
h_delta_flux_port_err = hdu_data.field(48)
h_gamma_flux_port = hdu_data.field(49)
h_gamma_flux_port_err = hdu_data.field(50)
oiii_4363_flux_port = hdu_data.field(51)
oiii_4363_flux_port_err = hdu_data.field(52)
heii_4685_flux_port = hdu_data.field(53)
heii_4685_flux_port_err = hdu_data.field(54)
ariv_4711_flux_port = hdu_data.field(55)
ariv_4711_flux_port_err = hdu_data.field(56)
h_beta_flux_port = hdu_data.field(57)
h_beta_flux_port_err = hdu_data.field(58)
oiii_4959_flux_port = hdu_data.field(59)
oiii_4959_flux_port_err = hdu_data.field(60)
oiii_5007_flux_port = hdu_data.field(61)
oiii_5007_flux_port_err = hdu_data.field(62)
hei_5876_flux_port = hdu_data.field(63)
hei_5876_flux_port_err = hdu_data.field(64)
oi_6300_flux_port = hdu_data.field(65)
oi_6300_flux_port_err = hdu_data.field(66)
nii_6548_flux_port = hdu_data.field(67)
nii_6548_flux_port_err = hdu_data.field(68)
h_alpha_flux_port = hdu_data.field(69)
h_alpha_flux_port_err = hdu_data.field(70)
nii_6584_flux_port = hdu_data.field(71)
nii_6584_flux_port_err = hdu_data.field(72)
sii_6717_flux_port = hdu_data.field(73)
sii_6717_flux_port_err = hdu_data.field(74)
sii_6731_flux_port = hdu_data.field(75)
sii_6731_flux_port_err = hdu_data.field(76)
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

heii_3203_flux_port_ext = np.zeros(len(EBV_excess))
oii_3726_flux_port_ext = np.zeros(len(EBV_excess))
oii_3729_flux_port_ext = np.zeros(len(EBV_excess))
neiii_3869_flux_port_ext = np.zeros(len(EBV_excess))
h_delta_flux_port_ext = np.zeros(len(EBV_excess))
h_gamma_flux_port_ext = np.zeros(len(EBV_excess))
oiii_4363_flux_port_ext = np.zeros(len(EBV_excess))
heii_4685_flux_port_ext = np.zeros(len(EBV_excess))
ariv_4711_flux_port_ext = np.zeros(len(EBV_excess))
h_beta_flux_port_ext = np.zeros(len(EBV_excess))
oiii_4959_flux_port_ext = np.zeros(len(EBV_excess))
oiii_5007_flux_port_ext = np.zeros(len(EBV_excess))
hei_5876_flux_port_ext = np.zeros(len(EBV_excess))
oi_6300_flux_port_ext = np.zeros(len(EBV_excess))
nii_6548_flux_port_ext = np.zeros(len(EBV_excess))
h_alpha_flux_port_ext  = np.zeros(len(EBV_excess))
nii_6584_flux_port_ext = np.zeros(len(EBV_excess))
sii_6717_flux_port_ext = np.zeros(len(EBV_excess))
sii_6731_flux_port_ext = np.zeros(len(EBV_excess))


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
    
    mw_A_heii_3203_flux_port_rest = mw_A_atmwlam[68]
    mw_A_oii_3726_flux_port_rest = mw_A_atmwlam[242]
    mw_A_oii_3729_flux_port_rest = mw_A_atmwlam[243]
    mw_A_neiii_3869_flux_port_rest = mw_A_atmwlam[290]
    mw_A_h_delta_flux_port_rest = mw_A_atmwlam[367]
    mw_A_h_gamma_flux_port_rest = mw_A_atmwlam[447]
    mw_A_oiii_4363_flux_port_rest = mw_A_atmwlam[454]
    mw_A_heii_4685_flux_port_rest = mw_A_atmwlam[562]
    mw_A_ariv_4711_flux_port_rest = mw_A_atmwlam[570]
    mw_A_h_beta_flux_port_rest= mw_A_atmwlam[620]
    mw_A_oiii_4959_flux_port_rest = mw_A_atmwlam[653]
    mw_A_oiii_5007_flux_port_rest = mw_A_atmwlam[669]
    mw_A_hei_5876_flux_port_rest = mw_A_atmwlam[959]
    mw_A_oi_6300_flux_port_rest = mw_A_atmwlam[1100]
    mw_A_nii_6548_flux_port_rest = mw_A_atmwlam[1183]
    mw_A_h_alpha_flux_port_rest  = mw_A_atmwlam[1187]
    mw_A_nii_6584_flux_port_rest = mw_A_atmwlam[1195]
    mw_A_sii_6717_flux_port_rest = mw_A_atmwlam[1239]
    mw_A_sii_6731_flux_port_rest = mw_A_atmwlam[1244]

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
    

    heii_3203_flux_port_deext_rest = 10.0**(mw_A_heii_3203_flux_port_rest/2.5) 
    oii_3726_flux_port_deext_rest = 10.0**(mw_A_oii_3726_flux_port_rest/2.5)
    oii_3729_flux_port_deext_rest = 10.0**(mw_A_oii_3729_flux_port_rest/2.5)
    neiii_3869_flux_port_deext_rest= 10.0**(mw_A_neiii_3869_flux_port_rest/2.5)
    h_delta_flux_port_deext_rest = 10.0**(mw_A_h_delta_flux_port_rest/2.5)
    h_gamma_flux_port_deext_rest = 10.0**(mw_A_h_gamma_flux_port_rest/2.5)
    oiii_4363_flux_port_deext_rest= 10.0**(mw_A_oiii_4363_flux_port_rest/2.5)
    heii_4685_flux_port_deext_rest = 10.0**(mw_A_heii_4685_flux_port_rest/2.5)
    ariv_4711_flux_port_deext_rest = 10.0**(mw_A_ariv_4711_flux_port_rest/2.5)
    h_beta_flux_port_deext_rest = 10.0**(mw_A_h_beta_flux_port_rest/2.5)
    oiii_4959_flux_port_deext_rest = 10.0**(mw_A_oiii_4959_flux_port_rest/2.5)
    oiii_5007_flux_port_deext_rest = 10.0**(mw_A_oiii_5007_flux_port_rest/2.5)
    hei_5876_flux_port_deext_rest = 10.0**(mw_A_hei_5876_flux_port_rest/2.5)
    oi_6300_flux_port_deext_rest = 10.0**(mw_A_oi_6300_flux_port_rest/2.5)
    nii_6548_flux_port_deext_rest = 10.0**(mw_A_nii_6548_flux_port_rest/2.5)
    h_alpha_flux_port_deext_rest  = 10.0**(mw_A_h_alpha_flux_port_rest/2.5)
    nii_6584_flux_port_deext_rest = 10.0**(mw_A_nii_6584_flux_port_rest/2.5)
    sii_6717_flux_port_deext_rest = 10.0**(mw_A_sii_6717_flux_port_rest/2.5)
    sii_6731_flux_port_deext_rest = 10.0**(mw_A_sii_6731_flux_port_rest/2.5)

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
    
    heii_3203_flux_port_ext[i] = heii_3203_flux_port[i]*heii_3203_flux_port_deext_rest
    oii_3726_flux_port_ext[i] = oii_3726_flux_port[i]*oii_3726_flux_port_deext_rest
    oii_3729_flux_port_ext[i] = oii_3729_flux_port[i]*oii_3729_flux_port_deext_rest
    neiii_3869_flux_port_ext[i] = neiii_3869_flux_port[i]*neiii_3869_flux_port_deext_rest
    h_delta_flux_port_ext[i] = h_delta_flux_port[i]*h_delta_flux_port_deext_rest
    h_gamma_flux_port_ext[i] = h_gamma_flux_port[i]*h_gamma_flux_port_deext_rest
    oiii_4363_flux_port_ext[i] = oiii_4363_flux_port[i]*oiii_4363_flux_port_deext_rest
    heii_4685_flux_port_ext[i] = heii_4685_flux_port[i]*heii_4685_flux_port_deext_rest
    ariv_4711_flux_port_ext[i] = ariv_4711_flux_port[i]*ariv_4711_flux_port_deext_rest
    h_beta_flux_port_ext[i] = h_beta_flux_port[i]*h_beta_flux_port_deext_rest
    oiii_4959_flux_port_ext[i] = oiii_4959_flux_port[i]*oiii_4959_flux_port_deext_rest
    oiii_5007_flux_port_ext[i] = oiii_5007_flux_port[i]*oiii_5007_flux_port_deext_rest
    hei_5876_flux_port_ext[i] = hei_5876_flux_port[i]*hei_5876_flux_port_deext_rest
    oi_6300_flux_port_ext[i] = oi_6300_flux_port[i]*oi_6300_flux_port_deext_rest
    nii_6548_flux_port_ext[i] = nii_6548_flux_port[i]*nii_6548_flux_port_deext_rest
    h_alpha_flux_port_ext[i]  = h_alpha_flux_port[i]*h_alpha_flux_port_deext_rest
    nii_6584_flux_port_ext[i] = nii_6584_flux_port[i]*nii_6584_flux_port_deext_rest
    sii_6717_flux_port_ext[i] = sii_6717_flux_port[i]*sii_6717_flux_port_deext_rest
    sii_6731_flux_port_ext[i] = sii_6731_flux_port[i]*sii_6731_flux_port_deext_rest

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

col36 = fits.Column(name = 'heii_3203_flux_port_ext', format = 'E', array=heii_3203_flux_port_ext)
col37 = fits.Column(name = 'heii_3203_flux_port_ext_err', format = 'E', array=heii_3203_flux_port_err)
col38 = fits.Column(name = 'oii_3726_flux_port_ext', format = 'E', array=oii_3726_flux_port_ext)
col39 = fits.Column(name = 'oii_3726_flux_port_ext_err', format = 'E', array=oii_3726_flux_port_err)
col40 = fits.Column(name = 'oii_3729_flux_port_ext', format = 'E', array=oii_3729_flux_port_ext)
col41 = fits.Column(name = 'oii_3729_flux_port_ext_err', format = 'E', array=oii_3729_flux_port_err)
col42 = fits.Column(name = 'neiii_3869_flux_port_ext', format = 'E', array=neiii_3869_flux_port_ext)
col43 = fits.Column(name = 'neiii_3869_flux_port_ext_err', format = 'E', array=neiii_3869_flux_port_err)
col44 = fits.Column(name = 'h_delta_flux_port_ext', format = 'E', array=h_delta_flux_port_ext)
col45 = fits.Column(name = 'h_delta_flux_port_ext_err', format = 'E', array=h_delta_flux_port_err)
col46 = fits.Column(name = 'h_gamma_flux_port_ext', format = 'E', array=h_gamma_flux_port_ext)
col47 = fits.Column(name = 'h_gamma_flux_port_ext_err', format = 'E', array=h_gamma_flux_port_err)
col48 = fits.Column(name = 'oiii_4363_flux_port_ext', format = 'E', array=oiii_4363_flux_port_ext)
col49 = fits.Column(name = 'oiii_4363_flux_port_ext_err', format = 'E', array=oiii_4363_flux_port_err)
col50 = fits.Column(name = 'heii_4685_flux_port_ext', format = 'E', array=heii_4685_flux_port_ext)
col51 = fits.Column(name = 'heii_4685_flux_port_ext_err', format = 'E', array=heii_4685_flux_port_err)
col52 = fits.Column(name = 'ariv_4711_flux_port_ext', format = 'E', array=ariv_4711_flux_port_ext)
col53 = fits.Column(name = 'ariv_4711_flux_port_ext_err', format = 'E', array=ariv_4711_flux_port_err)
col54 = fits.Column(name = 'h_beta_flux_port_ext', format = 'E', array=h_beta_flux_port_ext)
col55 = fits.Column(name = 'h_beta_flux_port_ext_err', format = 'E', array=h_beta_flux_port_err)
col56 = fits.Column(name = 'oiii_4959_flux_port_ext', format = 'E', array=oiii_4959_flux_port_ext)
col57 = fits.Column(name = 'oiii_4959_flux_port_ext_err', format = 'E', array=oiii_4959_flux_port_err)
col58 = fits.Column(name = 'oiii_5007_flux_port_ext', format = 'E', array=oiii_5007_flux_port_ext)
col59 = fits.Column(name = 'oiii_5007_flux_port_ext_err', format = 'E', array=oiii_5007_flux_port_err)
col60 = fits.Column(name = 'hei_5876_flux_port_ext', format = 'E', array=hei_5876_flux_port_ext)
col61 = fits.Column(name = 'hei_5876_flux_port_ext_err', format = 'E', array=hei_5876_flux_port_err)
col62 = fits.Column(name = 'oi_6300_flux_port_ext', format = 'E', array=oi_6300_flux_port_ext)
col63 = fits.Column(name = 'oi_6300_flux_port_ext_err', format = 'E', array=oi_6300_flux_port_err)
col64 = fits.Column(name = 'nii_6548_flux_port_ext', format = 'E', array=nii_6548_flux_port_ext)
col65 = fits.Column(name = 'nii_6548_flux_port_ext_err', format = 'E', array=nii_6548_flux_port_err)
col66 = fits.Column(name = 'h_alpha_flux_port_ext', format = 'E', array=h_alpha_flux_port_ext)
col67 = fits.Column(name = 'h_alpha_flux_port_ext_err', format = 'E', array=h_alpha_flux_port_err)
col68 = fits.Column(name = 'nii_6584_flux_port_ext', format = 'E', array=nii_6584_flux_port_ext)
col69 = fits.Column(name = 'nii_6584_flux_port_ext_err', format = 'E', array=nii_6584_flux_port_err)
col70 = fits.Column(name = 'sii_6717_flux_port_ext', format = 'E', array=sii_6717_flux_port_ext)
col71 = fits.Column(name = 'sii_6717_flux_port_ext_err', format = 'E', array=sii_6717_flux_port_err)
col72 = fits.Column(name = 'sii_6731_flux_port_ext', format = 'E', array=sii_6731_flux_port_ext)
col73 = fits.Column(name = 'sii_6731_flux_port_ext_err', format = 'E', array=sii_6731_flux_port_err)
col74 = fits.Column(name = 'balmer_decrement', format = 'E', array=balm_dec_obs)
col75 = fits.Column(name = 'balmer_decrement_ext', format = 'E', array=(h_alpha_flux_ext/h_beta_flux_ext))


#condense colummns, construct new table, write out results to binary fits
cols = fits.ColDefs([col1, col2,col3, col4, col5, col6, col7, col8, col9, col10, 
col11, col12, col13, col14, col15, col16, col17, col18, col19, col20, col21,col22,
col23,col24,col25,col26,col27,col28,col29,col30,col31,col32,col33,col34,col35, 
col36, col37, col38, col39, col40, col41, col42, col43, col44, col45, col46, col47,
col48, col49, col50, col51, col52, col53, col54, col55, col56, col57, col58, col59
, col60, col61, col62, col63, col64, col65, col66, col67, col68, col69, col70, col71, col72, col73, col74, col75])
#from astropy.table import Table
#tbhdu = fits.new_table(cols)
tbhdu = fits.BinTableHDU.from_columns(cols)
tbhdu.writeto('ECO_SDSS_full_dext.fits', clobber=True)
#tbhdu.writeto('RESOLVE_SDSS_dext_dwarf_SNR3.fits', clobber=True)
#####################################################################

