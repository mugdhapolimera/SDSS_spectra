# -*- coding: utf-8 -*-
"""
Created on Wed May  5 12:34:31 2021

@author: mugdhapolimera
"""

import pandas as pd
import numpy as np

agn = pd.read_csv("ECO+RESOLVE_AGN_list.csv")
agn.index = agn.name

opt = pd.read_csv("eco+resolve_emlineclass_dext_snr5_jhu.csv")
opt.index = opt.galname

jhu = pd.read_csv('ECO+RESOLVE_full_snr5_dext_jhu.csv')
jhu.index = jhu.name

resfull = pd.read_csv('C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO+RESOLVE_inobssample.csv')
resfull.index = resfull.name

#barro = pd.read_csv("Barro_inobssample.csv")
#barro.index = barro.name 

#Calculate OIII luminosity for optical emission line AGN
optagn = np.array(opt.galname[opt.defagn | opt.composite | opt.sftoagn | opt.agntosf])
o3flux = jhu.loc[optagn]['oiii_5007_flux'] * 1e-17 #ergs/s
vel = jhu.loc[optagn]['cz']
H0 = 70 #km/s/Mpc
Mpc = 3.086e24 #cm
d = (vel/H0) * Mpc #cm
o3lum = 4*3.14*(d**2) * o3flux

bol_lum = o3lum * 1e3
#agnbacklights = np.array(bol_lum.index[(bol_lum > 1e42) & \
#                                       (jhu.loc[bol_lum.index].logmstar > 9.5) &\
#                                       (barro.loc[bol_lum.index].MU_50 > 8.4)])

#dusty = ['rf0444', 'rs1036']
#agnbacklights = np.setdiff1d(agnbacklights, dusty)

#jhu.loc[agnbacklights][['radeg','dedeg']].to_csv("RESOLVE_AGN_backlights.csv")


