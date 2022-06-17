# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 23:08:09 2021

@author: mugdhapolimera
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams.update({'axes.linewidth': 2})
matplotlib.rcParams.update({'lines.linewidth': 2})
import scipy
from scipy.io.idl import readsav
from collections import OrderedDict
from matplotlib import cm
import math
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.cbook import get_sample_data
import matplotlib.image as mpimg
import matplotlib.colors as mpcolors

def density_estimation(m1, m2):
    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]                                                     
    positions = np.vstack([X.ravel(), Y.ravel()])                                                       
    values = np.vstack([m1, m2])                                                                        
    kernel = scipy.stats.gaussian_kde(values)                                                                 
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z 

resolve = pd.read_csv("RESOLVE_inobssample.csv")
resolve.index = resolve.name
internal = readsav("resolvecatalog.dat")
internalphot = readsav("resolvecatalogphot.dat")
res_el = pd.read_csv('RESOLVE_hasnr5_inobssamplejhu.csv')
res_el.index = res_el.name
res_bpt = pd.read_csv('RESOLVE_bpt1snr5_inobssamplejhu.csv')
res_bpt.index = res_bpt.name

names, resndx, catndx = np.intersect1d(resolve.name, internal.name, return_indices = True)

blueflag = internal.inobssample & ((internal.blue) | (internal.gemdone) | \
        (internal.koalablue) | (internal.saltlsblue))
redflag = internal.inobssample & ((internal.red) | (internal.gemdone) | \
        (internal.koalared) | (internal.saltls))
broadflag = internal.inobssample & (internal.broad)
u_r = resolve['modelu_rcorr']
xmin = 7.5
xmax = 11.5
ymin = 0
ymax = 3

X,Y,Z = density_estimation(resolve.logmstar,u_r)
fig,ax = plt.subplots()#(figsize=(8,8))
Z = Z/Z.max()
lvls = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
def truncate_colormap(cmap, minval=0, maxval=0.75, n=150):
  	new_cmap = mpcolors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
  	return new_cmap
sf_colors_map = truncate_colormap(cm.gray,n=11)
nbins = 20
plt.contour(X, Y, Z, levels = lvls, cmap=sf_colors_map, zorder = 1)
#plt.imshow(np.rot90(Z), cmap='bone_r',                                                    
#          extent=[xmin, xmax, ymin, ymax], interpolation='gaussian')
#plt.clim(0,1.8)
#plt.contour(X, Y, Z, cmap='summer')

broad = np.intersect1d(names, internal.name[np.where(internal.broad)])
blue = np.intersect1d(names, internal.name[np.where(blueflag)])
red = np.intersect1d(names, internal.name[np.where(redflag)])

newspec = list(broad) + list(blue) + list(red)
newspec = np.unique(newspec)
plt.plot(resolve.logmstar.loc[newspec], 
        u_r.loc[newspec], 'k*', markersize = 8, label = 'New RESOLVE Spectra'),
plt.plot(res_el.logmstar.loc[newspec], 
        res_el.loc[newspec].modelu_rcorr, 'o', mfc = 'none',
        color = 'lime', markersize = 10, label = 'Expected Emission-Line Sample'),
#plt.plot(res_bpt.logmstar.loc[newspec], 
#        res_bpt.loc[newspec].modelu_rcorr, 'k*', markersize = 8, label = 'SDSS BPT Emission-Line Sample'),
plt.legend()
plt.xlabel('log(Stellar Mass/M$_\odot$)')
plt.ylabel('(u-r) colour')
#agn = pd.read_csv("ECO+RESOLVE_AGN_list.csv")
#agn.index = agn.name
#agn = agn[(agn.agntype != 'xrayagn') & (agn.agntype != 'midiragn')] #agn[agn.bptagn | agn.agntosf | agn.bptcomposite | agn.sfingagn]
#sfagn = np.intersect1d(agn.name[agn.agntype == 'sfingagn'], resolve.name)
#bptagn = np.intersect1d(agn.name[agn.agntype == 'bptagn'], resolve.name)
#compositeagn = np.intersect1d(agn.name[agn.agntype == 'bptcomposite'], resolve.name)
#
#plt.plot(resolve.logmstar.loc[sfagn], 
#        u_r.loc[sfagn], 'bs', markersize = 8)#, label = 'SFing-AGN'),
#plt.plot(resolve.logmstar.loc[bptagn], 
#        u_r.loc[bptagn], 'rs', markersize = 8)#, label = 'SFing-AGN'),
#plt.plot(resolve.logmstar.loc[compositeagn], 
#        u_r.loc[compositeagn], 'ms', markersize = 8)#, label = 'SFing-AGN'),
#
#
#midiragn = pd.read_csv(r'C:\Users\mugdhapolimera\github\SDSS_spectra\mid_ir\RESOLVE_WISE_AGN.csv')
#midiragn.index = midiragn.name
#iragn = np.intersect1d(midiragn.name, resolve.name)
#plt.plot(resolve.logmstar.loc[iragn], 
#        u_r.loc[iragn], 'p', color = 'orange', markersize = 8)#, label = 'SFing-AGN'),
#
