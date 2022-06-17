# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 11:02:49 2020

@author: mugdhapolimera
https://joseph-long.com/writing/from-sky-coordinates-to-pixels-and-back/
"""

import numpy as np
import pandas as pd
from scipy.io import readsav
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import os
#from astropy.wcs import WCS as wcs 
from astropy import wcs
from astropy.coordinates import SkyCoord as sky
from astropy import units as u
from matplotlib.collections import PatchCollection

#os.chdir(r'F:\mugdhapolimera\Download\/')
image = r'C:\Users\mugdhapolimera\github\xray\catalog_matching\radio\postage (4).fits'#_NVSS_.fits'
#cutout_135.2130_1.1612.fits'
name = image#'rs0010.jfif'
#jpg = plt.imread(name)
jpg = fits.open(image)[0].data
if jpg.ndim > 3:
    jpg = jpg[0][0]

#jpg2 = np.flipud(jpg)
cubehdu = fits.open(image)[0].header
cube = cubehdu#fits.open(image)[0].data#.flatten()

w = wcs.WCS(cubehdu)
w = w.dropaxis(2)
w = w.dropaxis(2)
xy = np.meshgrid(np.arange(50), np.arange(50))
coords = w.all_pix2world(xy[0], xy[1],0)
fig = plt.figure()
ax = plt.subplot(projection=w)

from astropy.visualization import make_lupton_rgb
#rgb_default = make_lupton_rgb(cube[2], cube[1], cube[0], Q = 10, stretch = 0.5)
ax.imshow(jpg/np.max(jpg), origin='lower')#, cmap=plt.cm.viridis)
#ax.scatter(sdsspix[0][0], sdsspix[0][1],color = 'm')
plt.xlabel('RA', fontsize = 22)
plt.ylabel('Dec', fontsize = 22)
#plt.ylim(firstlims)
if image == r'C:\Users\mugdhapolimera\github\xray\catalog_matching\radio\rs1038':
    
    yfirstlims = plt.gca().get_ylim()
    ywcslims = w.all_pix2world(yfirstlims[0],yfirstlims[1],0)
    xfirstlims = plt.gca().get_xlim()
    xwcslims = w.all_pix2world(xfirstlims[0],xfirstlims[1],0)
    
else:
    (xmin, xmax) = w.all_world2pix(xwcslims[0], xwcslims[1],0)
    (ymin, ymax) = w.all_world2pix(ywcslims[0], ywcslims[1],0)
    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)