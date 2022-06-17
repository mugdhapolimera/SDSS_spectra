# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 12:26:01 2020

@author: mugdhapolimera
"""

import pandas as pd
import numpy
import numpy as np
from scipy.io import readsav
rescat = readsav('resolvecatalog.dat')
resphotcat = readsav('resolvecatalogphot.dat')
jhuflag = pd.read_csv('resolve_emlineclass_dext_snr5_jhu.csv')
jhuflag.index = jhuflag.galname
sfingagn = np.array(jhuflag.galname[jhuflag.sftoagn])
sfingagnndx = [x for x in range(len(rescat.name)) if rescat.name[x] in sfingagn]
print rescat.ifusb[sfingagnndx]
print np.median(rescat.ifusb[sfingagnndx])
print np.median(rescat.vhel[sfingagnndx]/3e5)
print (rescat.vhel[sfingagnndx]/3e5)

jhu = pd.read_csv("RESOLVE_full_snr5_dext_jhu.csv")
jhu.index = jhu.name
print jhu.loc[sfingagn].cz/3e5
print jhu.loc[sfingagn][['radeg','dedeg']]
import astropy.io.fits
from astropy.io import fits
spec= fits.open("F:\mugdhapolimera\Downloads\spec-0470-51929-0336.fits")
import matplotlib.pyplot as plt

plt.figure()
plt.plot(10**spec[1].data.loglam, spec[1].data.flux)
plt.plot(10**spec[1].data.loglam, spec[1].data.model)

lam = 10**spec[1].data.loglam
flux = spec[1].data.model

lam = lam/10
print rescat.ifusb[np.where(rescat.name == 'rs0010')]
flux = flux*1e-17
flux = flux*10
np.savetxt(r"C:\Users\mugdhapolimera\github\SDSS_spectra\rs0010spec.txt",zip(lam,flux))
#rescat.r50[np.where(rescat.name == 'rs0010')]
#resphot.r50[np.where(rescat.name == 'rs0010')]
#resphotcat.r50[np.where(rescat.name == 'rs0010')]
#resphotcat.reff[np.where(rescat.name == 'rs0010')]
#resphotcat.re[np.where(rescat.name == 'rs0010')]
#rescat.re[np.where(rescat.name == 'rs0010')]
#rescat.reff[np.where(rescat.name == 'rs0010')]
#resphotcat.radr50p[np.where(rescat.name == 'rs0010')]
#resphotcat.b_a[np.where(rescat.name == 'rs0010')]