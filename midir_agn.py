# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 14:54:24 2019

@author: mugdhapolimera

Mid-IR AGN Selection using WISE IR colours

as prescribed by Sartori et al. 2015 in the method of Jarrett et al 2015
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from scipy.io.idl import readsav
#os.chdir('C:\Users\mugdhapolimera\github\SDSS_spectra')
def stern(x):
    return 0.8*np.ones(len(x))

def jarretx(y):
    return [2.2*np.ones(len(y)), 4.2*np.ones(len(y))]

def jarrety(x):
    return [1.7*np.ones(len(x)), 0.1*x+0.38]
os.chdir('/afs/cas.unc.edu/users/m/u/mugpol/github/SDSS_spectra/')
data = pd.read_csv('ECO+RESOLVE_filter_new.csv')

os.chdir('/afs/cas.unc.edu/users/m/u/mugpol/github/izi/')
catalog = readsav('resolvecatalogphot.dat')

catalog2 = readsav('resolvecatalog.dat')

w12 = catalog['mw1'] - catalog['mw2']
w23 = catalog['mw2'] - catalog['mw3']

plt.figure()
plt.plot(w23,w12,'o')
xaxis = np.linspace(min(w23),max(w23))
#yaxis = np.linspace(min(w12), max(w12))
yaxis = np.linspace(jarrety(np.array([2.2]))[1],1.7)
plt.plot(xaxis, stern(xaxis), 'k--')

xaxis = np.linspace(2.2,4.2)
plt.plot(jarretx(yaxis)[0], yaxis, 'k', jarretx(yaxis)[1], yaxis, 'k')
plt.plot(xaxis, jarrety(xaxis)[0], 'k', xaxis, jarrety(xaxis)[1],'k')

midiragn = (w12 >= 0.8) | ((w23 > 2.2) & (w23 < 4.2) & (w12 < 1.7) & 
                            (0.1*w23 + 0.38 < w12))
plt.plot(w23[midiragn],w12[midiragn],'rs')

for i in range(len(catalog2.name[midiragn])):
    print "'"+catalog2.name[midiragn][i]+"' ,",