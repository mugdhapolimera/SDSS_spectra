# -*- coding: utf-8 -*-
'''
Created on Tue Apr 21 10:56:12 2020

@author: mugdhapolimera

Plot ECO galaxies with official aspect ratio
'''
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.coordinates import SkyCoord
import matplotlib as mpl
import astropy.units as u

fulleco = pd.read_csv('ECO_live22Oct2018.csv')
voleco = pd.read_csv('ECO_inobssample.csv')
eleco = pd.read_csv('ECO_snr5_master_hasnr5.csv')
seleco = pd.read_csv('ECO_snr5_master_bary.csv')


fig = plt.figure()
ax = fig.add_subplot(111)

ras = fulleco.radeg
decs = fulleco.dedeg
c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree)
ras = c.ra.hour

#ax.plot(ras, decs, 'ko', ms=3, mfc = 'none', label='Galaxy')
ax.set_xlabel('RA (hours)', fontsize = 15)
ax.set_ylabel('Dec (degrees)', fontsize = 15)
ax.set_aspect(aspect=0.0667)
#ax.set_title('Galaxies in ECO')
#ax.hlines(5., min(ras), max(ras), colors='red', linestyles='dashed', lw=3)
#cmap = plt.get_cmap('rainbow',100)
#norm = mpl.colors.Normalize()
#sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
#sm.set_array([])
#cbar = plt.colorbar(sm)
#cbar.set_label('Normalized Density of Galaxies', rotation=270, labelpad = 35)
plt.xlim(min(ras), max(ras))
plt.ylim(min(decs), 50)#max(decs))

ras = voleco.radeg
decs = voleco.dedeg
c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree)
ras = c.ra.hour

ax.scatter(ras, decs, s=3, color = 'gray', label='Mass-limited Sample')

ras = eleco.radeg
decs = eleco.dedeg
c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree)
ras = c.ra.hour

ax.scatter(ras, decs, s=3, color = 'orange', label='Emission Line Sample')

ras = seleco.radeg
decs = seleco.dedeg
c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree)
ras = c.ra.hour

ax.scatter(ras, decs, s=3, color = 'black', label='Strong Emission Line Sample')

#fontP = FontProperties()
#fontP.set_size('small')
lgnd = plt.legend(bbox_to_anchor=(1.01,1), loc="upper left", fontsize = 15,
                  title = 'ECO')
for handle in lgnd.legendHandles:
    handle.set_sizes([20.0])
lgnd.get_title().set_fontsize('15')
#plt.setp(plt.gca().get_legend().get_texts(), fontsize='15') #legend 'list' fontsize
plt.show()
