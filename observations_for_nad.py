# -*- coding: utf-8 -*-
"""
Created on Wed May  5 20:13:09 2021

@author: mugdhapolimera
"""

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import os
from scipy.io.idl import readsav
import numpy as np

os.chdir('C:\Users\mugdhapolimera\github\SDSS_spectra/')
resolve = pd.read_csv('RESOLVE_full_raw.csv')
resolve.index = resolve.name
internal = readsav('resolvecatalog.dat')
internalphot = readsav('resolvecatalogphot.dat')
barro = pd.read_csv("Barro_inobssample.csv")
barro.index = barro.name 

agn = pd.read_csv("RESOLVE_AGN_list.csv")
agn.index = agn.name

nad_targets_a = internal.name[(internal.broad == 1) & (internal.obstag == 'A')]
nad_targets_b = internal.name[(internal.broad == 1) & (internal.obstag == 'B')]

keys = ['sfingagn', 'bptagn', 'bptcomposite', 'midiragn', 'xrayagn', 'agntosf'] 
    
marker = {'agntosf': '^', 'bptcomposite': 's', 'bptagn': 's', 
          'sfingagn': 's', 'xrayagn': 's' ,'midiragn': 'v'}
colors = {'agntosf': 'c', 'bptcomposite': 'brown', 'bptagn': 'brown', 
          'sfingagn': 'b', 'xrayagn': 'k', 'midiragn': 'lime'}
labels = {'agntosf': 'Low-SII AGN', 'bptcomposite': 'BPT AGN', 'bptagn': 'BPT AGN', 
          'sfingagn': 'SF-AGN', 'xrayagn': 'X-ray AGN', 'midiragn': 'Mid-IR AGN'}
zo = 0
plt.figure()
plt.plot(internalphot.absmagr, internal.ifusb, 'o', color = 'gray', zorder = zo)

for key in keys:
    zo += 1
    if key != 'bptagn':
        sel = [x for x in range(len(internal.name)) if internal.name[x] in list(agn.name[agn[key]])]
        if len(sel) > 0:
                plt.scatter(internalphot.absmagr[sel],internal.ifusb[sel], 
                     marker = marker[key], s = 150,
                     color = colors[key], label = labels[key], zorder = zo)
plt.legend(fontsize = 15)
key = 'bptagn'
zo +=1
sel = [x for x in range(len(internal.name)) if internal.name[x] in list(agn.name[agn[key]])]
if len(sel) > 0:
        plt.scatter(internalphot.absmagr[sel],internal.ifusb[sel], 
             marker = marker[key], s = 200,
             color = colors[key], label = labels[key], zorder = 1)
plt.xlim(-15,-24)
plt.ylim(24,16)
plt.xlabel(r'Absolute r-band mag')
plt.ylabel(r'Central Surface Brightness (mag/arcsec$^2$)')
#check running logs and see if the data for these are fine
#maybe exposure time can hint at s/n?

zo = 0
plt.figure()
plt.plot(internalphot.absmagr, barro.loc[internal.name].MU_50, 'o', color = 'gray', zorder = zo)

for key in keys:
    zo += 1
    if key != 'bptagn':
        sel = [x for x in range(len(internal.name)) if internal.name[x] in list(agn.name[agn[key]])]
        if len(sel) > 0:
                plt.scatter(internalphot.absmagr[sel],barro.loc[internal.name[sel]].MU_50, 
                     marker = marker[key], s = 150,
                     color = colors[key], label = labels[key], zorder = zo)
plt.legend(fontsize = 15)
key = 'bptagn'
zo +=1
sel = [x for x in range(len(internal.name)) if internal.name[x] in list(agn.name[agn[key]])]
if len(sel) > 0:
        plt.scatter(internalphot.absmagr[sel],barro.loc[internal.name[sel]].MU_50, 
             marker = marker[key], s = 200,
             color = colors[key], label = labels[key], zorder = 1)
plt.xlim(-15,-24)
plt.ylim(5,10)
plt.xlabel(r'Absolute r-band mag', fontsize = 20)
plt.ylabel(r'Compactness metric $\Delta\mu_{50}$ or $\Delta\Sigma_e$', fontsize = 20)
