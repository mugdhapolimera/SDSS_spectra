# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 11:02:49 2020

@author: mugdhapolimera
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
image = 'cutout_135.2130_1.1612.fits'
name = 'rs0010.jfif'
jpg = plt.imread(name)
jpg2 = np.flipud(jpg)
cubehdu = fits.open(image)[0].header
cube = fits.open(image)[0].data#.flatten()

w = wcs.WCS(cubehdu)
w = w.dropaxis(2)
xy = np.meshgrid(np.arange(50), np.arange(50))
coords = w.all_pix2world(xy[0], xy[1],0)
fig = plt.figure()
ax = plt.subplot(projection=w)

from astropy.visualization import make_lupton_rgb
rgb_default = make_lupton_rgb(cube[2], cube[1], cube[0], Q = 10, stretch = 0.5)
ax.imshow(jpg2, origin='lower')#, cmap=plt.cm.viridis)
#ax.scatter(sdsspix[0][0], sdsspix[0][1],color = 'm')
plt.xlabel('RA', fontsize = 22)
plt.ylabel('Dec', fontsize = 22)
