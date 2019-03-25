# -*- coding: utf-8 -*-
"""
Created on Thu Mar 07 16:32:49 2019

@author: mugdhapolimera
"""

import numpy as np
import pandas as pd
import astropy
from astropy.io import fits
from astropy import units as u
from matplotlib import pyplot as plt
from astropy.visualization import quantity_support
import os 
from os.path import join
from specutils.fitting import fit_generic_continuum
from specutils import Spectrum1D
from specutils import SpectralRegion
from specutils.analysis import equivalent_width
quantity_support()

from pysynphot import observation
from pysynphot import spectrum
 
def rebin_spec(wave, specin, wavnew):
    spec = spectrum.ArraySourceSpectrum(wave=wave, flux=specin)
    f = np.ones(len(wave))
    filt = spectrum.ArraySpectralElement(wave, f, waveunits='angstrom')
    obs = observation.Observation(spec, filt, binset=wavnew)
 
    return obs.binflux

path = r'F:\mugdhapolimera\Documents\UNC\Research\Data\RESOLVE\SDSS_spectra_fits'
os.chdir(path)

filenames = []
for (dirpath, dirnames, filelist) in os.walk(path):
    filenames += [join(dirpath, file) for file in filelist]
speclistname = r'F:\mugdhapolimera\Documents\UNC\Research\Data\RESOLVE\SDSS_spectra_fits\speclist.txt'
speclist = np.loadtxt(speclistname,dtype = str)
resname = pd.read_csv(filenames[0])
datalist = pd.Series(index = list(resname['name']))#,columns=['filename'])

for i in range(len(speclist)): 
    datalist.loc[resname['name'][i]] = (join(path,speclist[i]).replace("\\","/")).replace("/","\\")
    
resolve = pd.read_pickle(r'C:\Users\mugdhapolimera\github\SDSS_spectra\RESOLVE_filter.pkl')
filenames = []
#for name in resolve['NAME']:
 #   filenames.append(datalist[name])

for gal in ['rf0031']:#datalist.index:
    if gal in resolve.index:
        print (gal)
        filename = datalist[gal] #'F:\\mugdhapolimera\\Documents\\UNC\\Research\\Data\\RESOLVE\\SDSS_spectra_fits\\0390\\spec-0390-51900-0300.fits'#
        # The spectrum is in the second HDU of this file.
        f = fits.open(filename)
        specdata = f[1].data # doctest: +REMOTE_DATA
        f.close()
    
        lamb_full = 10**specdata['loglam'] * u.AA # doctest: +REMOTE_DATA
        lamb = lamb_full[lamb_full < 6900*u.AA]
        
        flux = specdata['flux'] * 10**-17 * u.Unit('erg cm-2 s-1 AA-1') # doctest: +REMOTE_DATA
        flux = flux[lamb_full < 6900*u.AA]
        spec = Spectrum1D(spectral_axis = lamb, flux = flux) # doctest: +REMOTE_DATA
        cont_norm_spec = spec / fit_generic_continuum(spec)(spec.spectral_axis) # doctest: +REMOTE_DATA
        
        wave = np.array(cont_norm_spec.wavelength)
        flux = np.array(cont_norm_spec.flux)
        new_wave = wave/(1+(resolve.loc[gal].cz/3e5))
        new_flux = flux*(1+(resolve.loc[gal].cz/3e5))#rebin_spec(new_wave,flux,wave)
        new_spec = Spectrum1D(spectral_axis = new_wave* u.AA , flux = new_flux *u.Unit('1')) # doctest: +REMOTE_DATA
        plt.figure()
        plt.plot( new_wave,new_flux,'b')
        
        center = 6562#4685
        lower = center - 10
        upper =  center + 10
        plt.figure()
        plt.plot(new_wave,new_flux,'b')
        plt.xlim(lower,upper)
        print(equivalent_width(new_spec, regions=SpectralRegion(lower*u.AA, upper*u.AA)))
        
        center = 4711
        lower = center - 10
        upper =  center + 10
        plt.figure()
        plt.plot(new_wave,new_flux,'b')
        plt.xlim(lower,upper)
        print(equivalent_width(new_spec, regions=SpectralRegion(lower*u.AA, upper*u.AA)))
        
        #from spectres import spectres
        #x = spectres(wave, new_wave, flux)
        #print (x)