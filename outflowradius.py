# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 19:18:05 2021

@author: mugdhapolimera
"""

import numpy as np
import pandas as pd

inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO+RESOLVE_snr5_dext_jhu.csv'
df = pd.read_csv(inputfile)
df.index = df.name
flags = pd.read_csv(r'ECO\SEL\eco+resolve_emlineclass_dext_snr5_jhu.csv')
flags.index= flags.galname

#[O III] fluxes of SF-AGN; units = erg /s/cm^2
o3 = df['oiii_5007_flux'][flags.sftoagn] * 10**-17

H0 = 70 #km/s/Mpc
Mpc = 3.086e24 #Mpc to cm conversion

#distance to galaxy in cm
d = df['cz'][flags.sftoagn]/H0 * Mpc

#[O III] luminosity in erg/s
L_o3 = 4*3.14*(d**2)*o3

#NLR extent: log(R_NLR/1000 pc)
#logRkpc = (0.22* np.log(L_o3/1e42))+ 3.76 #Greene11
logR = (0.25*np.log10((L_o3/1e42))) + 3.746 #Liu13

Rkpc = (10**logR)/1000
#Rkpc = (np.exp(logRkpc))
print Rkpc