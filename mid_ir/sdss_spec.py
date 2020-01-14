# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 15:29:14 2019

@author: mugdhapolimera

Plotting SDSS spectra of 3 mid-IR AGN nuggets
"""

import astropy.io.fits as fits
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np

rs0814 = Table.read('spec-2880-54509-0230.fits')
plt.figure()
plt.plot(10**rs0814['loglam'], rs0814['flux'])

rf0073 = Table.read('spec-0395-51783-0525.fits')
plt.figure()
plt.plot(10**rf0073['loglam'], rf0073['flux'])

rf0006 = Table.read('spec-0390-51900-0282.fits')
plt.figure()
plt.plot(10**rf0006['loglam'], rf0006['flux'])

def func(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))    
def composite_spectrum(x, # data
                       a1, x01, sigma1, # 1st line
                       a2, x02, sigma2, # 2nd line
                       a3, x03, sigma3, # 3rd line
                       a4, x04, sigma4): # 4th line    
    return (func(x, a1, x01, sigma1) \
            + func(x, a2, x02, sigma2) \
            + func(x, a3, x03, sigma3)\
            + func(x, a4, x04, sigma4))
#                    + func(x, a3, x03, sigma3))
wavelengths = {'nii6548' : 6548.050, 
               'halpha' : 6562.819, 'nii6584' : 6583.46}

newlam = 10**rs0814['loglam']/(1+0.018)
totalflux = rs0814['flux']*(1+0.018)
from scipy.optimize import curve_fit
x2 = newlam[(newlam>6540) & (newlam<6590)]
y2 = totalflux[(newlam>6540) & (newlam<6590)]/1e-16
y2 = y2 - np.min(y2)
y2 = y2/1e19
#sig_narrow = 1.69/2.354 #GEMINI
sig_narrow = 69/3e5*6563 #SDSS

guess = [1000, 6563, sig_narrow, 500, 6584, sig_narrow,
         200, 6548, sig_narrow, 150, 6563, 4]

popt, pcov = curve_fit(composite_spectrum, x2, y2, p0 = guess)
plt.figure()
plt.plot(x2,y2)
plt.plot(x2, composite_spectrum(x2, *popt), 'k', label='Total fit')
plt.plot(x2, func(x2, *popt[:3]), c='g', lw = 3,
         label='Narrow Halpha')
plt.plot(x2, func(x2, *popt[-3:]), c='r', 
         label='Broad Halpha')
FWHM = round(2*np.sqrt(2*np.log(2))*popt[-1],4)
#plt.axvspan(popt[-2]-FWHM/2, popt[-2]+FWHM/2, 
#            facecolor='g', alpha=0.3, label='FWHM = %s'%(FWHM))
plt.legend(fontsize=10)
#plt.ylim(0,1.1)
