# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 14:54:24 2019

@author: mugdhapolimera

Mid-IR AGN Selection using WISE IR colours

as prescribed by Sartori et al. 2015 using the method of Jarrett et al 2015, 
Stern et al 2012; and the method of Satyapal et al. 2018 
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.io.idl import readsav
import sys
#os.chdir('C:\Users\mugdhapolimera\github\SDSS_spectra')
def stern(x):
    return 0.8*np.ones(len(x))

def jarretx(y):
    return [2.2*np.ones(len(y)), 4.2*np.ones(len(y))]

def jarrety(x):
    return [1.7*np.ones(len(x)), 0.1*x+0.38]

def satyapalx(x):
    return 0.52*np.ones(len(x))

def satyapaly(x):
    return 5.78*x -24.50
if sys.platform == 'linux2':
    os.chdir('/afs/cas.unc.edu/users/m/u/mugpol/github')
else:
    os.chdir('C:/Users/mugdhapolimera/github')
data = pd.read_csv('SDSS_spectra/ECO+RESOLVE_filter_new.csv')
catalog = readsav('izi/resolvecatalogphot.dat')

catalog2 = readsav('izi/resolvecatalog.dat')

baderr = np.isnan(catalog.emw1) | np.isnan(catalog.emw2) | np.isnan(catalog.emw3) | \
        np.isnan(catalog.emw4)
badphot = (catalog.mw1 == 0.0) & (catalog.mw2 == 0.0) & (catalog.mw3 == 0.0) & \
            (catalog.mw4 == 0.0)
bad = baderr | badphot
#catalog = catalog[~bad]
w12 = catalog['mw1'][~bad] - catalog['mw2'][~bad]
w23 = catalog['mw2'][~bad] - catalog['mw3'][~bad]
w12_err = np.sqrt(catalog['emw1'][~bad]**2 + catalog['emw2'][~bad]**2)
w23_err = np.sqrt(catalog['emw2'][~bad]**2 + catalog['emw3'][~bad]**2)
resname = catalog2['name'][~bad]
plt.figure()
plt.plot(w23,w12,'o')
xaxis = np.linspace(min(w23),max(w23))
#yaxis = np.linspace(min(w12), max(w12))
yaxis = np.linspace(jarrety(np.array([2.2]))[1],1.7)
plt.plot(xaxis, stern(xaxis), 'k--')
plt.plot(xaxis, satyapalx(xaxis), 'k-.')
xaxis = np.linspace(4.3,max(w23))
plt.plot(xaxis, satyapaly(xaxis), 'k-.')

xaxis = np.linspace(2.2,4.2)
plt.plot(jarretx(yaxis)[0], yaxis, 'k', jarretx(yaxis)[1], yaxis, 'k')
plt.plot(xaxis, jarrety(xaxis)[0], 'k', xaxis, jarrety(xaxis)[1],'k')
plt.xlabel('W2 - W3')
plt.ylabel('W1 - W2')
#midiragn = ((w12 >= 0.8) | ((w23 > 2.2) & (w23 < 4.2) & (w12 < 1.7) & 
#                            (0.1*w23 + 0.38 < w12)) | 
#           (w12 >= 0.52) & (w12 >= (5.78*w23) -24.50))

midiragn = ((w12-w12_err >= 0.8) | ((w23 > 2.2) & (w23 < 4.2) & \
             (w12-w12_err < 1.7) & (0.1*w23 + 0.38 < w12-w12_err)) | \
           (w12-w12_err >= 0.52) & (w12-w12_err >= (5.78*w23) -24.50))

#plt.plot(w23[midiragn],w12[midiragn],'rs')
plt.errorbar(w23[midiragn],w12[midiragn],fmt = 'rs', xerr = w23_err[midiragn],
             yerr = w12_err[midiragn])
#plt.ylim(-6.0,10)
#To print names of mid-IR AGN
for i in range(len(catalog2.name[~bad][midiragn])):
    print "'"+catalog2.name[~bad][midiragn][i]+"' ,",


#Flags to check SFing-AGN    
flags = pd.read_csv('SDSS_spectra/resolve_emlineclass_full_snr5.csv')
flags.index = flags.galname
    
sfagn = list(flags.galname.iloc[np.where(flags.sftoagn)])
sfagnndx = [x for x in range(len(resname)) if resname[x] in sfagn]
plt.plot(w23[sfagnndx],w12[sfagnndx],'gs')

sfagnmidir = [resname[midiragn][x] for x in range(len(resname[midiragn])) if resname[midiragn][x] in sfagn]
print sfagnmidir
