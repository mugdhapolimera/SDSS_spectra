from __future__ import print_function

#import pyfits
from astropy.table import Table
from scipy import ndimage
from scipy.io.idl import readsav
from scipy.optimize import curve_fit
import numpy as np
from time import clock
import glob


from numpy import pi
from numpy.ma import median
from matplotlib import pyplot as plt
import sys
import os
r = readsav('/srv/one/resolve/database_internal/merged_idl_catalog/stable/resolvecatalog.dat')
sel= np.where((r.inobssample == 1))
#dwarf
#sel= np.where((r.inobssample == 1) & (np.log10(r.mstars) < 9.5))
galname= r.name[sel]
ra = r.ra[sel]
dec = r.dec[sel]

create = np.array([galname, ra, dec])
#oschdir('/afs/cas.unc.edu/users/m/u/mugpol/Desktop/')
np.savetxt('fullRESOLVEcoords.txt', np.transpose(create), fmt ="%s")
