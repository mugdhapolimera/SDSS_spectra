# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 06:17:03 2022

@author: mugdhapolimera
"""

import pandas as pd
import numpy as np
from scipy.io.idl import readsav

#Baryonic-mass limited sample of RESOLVE and ECO
fullres = pd.read_csv("RESOLVE_barysample.csv")
dwarffullres = fullres.loc[fullres.logmstar < 9.5]
fulleco = pd.read_csv("ECO_barysample.csv")
dwarffulleco = fulleco.loc[fulleco.logmstar < 9.5]
print('Number of Galaxies in Baryonic mass-limited sample:')
print('RESOLVE :', len(fullres))
print('RESOLVE dwarfs:', len(dwarffullres))
print('Number of Galaxies:')
print('ECO :', len(fulleco))
print('ECO dwarfs:', len(dwarffulleco))

###############################################################################
#SEL sample of RESOLVE and ECO from MPA-JHU
selres = pd.read_csv("RESOLVE_full_snr5_dext_jhu.csv")
dwarfselres = selres.loc[selres.logmstar < 9.5]
seleco = pd.read_csv("ECO/SEL/ECO_full_snr5_dext_jhu.csv")
dwarfseleco = seleco.loc[seleco.logmstar < 9.5]
print('Number of Galaxies in JHU SEL sample:')
print('RESOLVE :', len(selres))
print('RESOLVE dwarfs:', len(dwarfselres))
print('Number of Galaxies:')
print('ECO :', len(seleco))
print('ECO dwarfs:', len(dwarfseleco))

#SEL sample of RESOLVE and ECO from Portsmouth
selres = pd.read_csv("RESOLVE_full_snr5_dext_port.csv")
dwarfselres = selres.loc[selres.logmstar < 9.5]
seleco = pd.read_csv("ECO/SEL/ECO_full_snr5_dext_port.csv")
dwarfseleco = seleco.loc[seleco.logmstar < 9.5]
print('Number of Galaxies in Portsmouth SEL sample:')
print('RESOLVE :', len(selres))
print('RESOLVE dwarfs:', len(dwarfselres))
print('Number of Galaxies:')
print('ECO :', len(seleco))
print('ECO dwarfs:', len(dwarfseleco))

#SEL sample of RESOLVE and ECO from NSA
selres = pd.read_csv("RESOLVE_full_snr5_dext_nsa.csv")
dwarfselres = selres.loc[selres.logmstar < 9.5]
seleco = pd.read_csv("ECO/SEL/ECO_full_snr5_dext_nsa.csv")
dwarfseleco = seleco.loc[seleco.logmstar < 9.5]
print('Number of Galaxies in NSA SEL sample:')
print('RESOLVE :', len(selres))
print('RESOLVE dwarfs:', len(dwarfselres))
print('Number of Galaxies:')
print('ECO :', len(seleco))
print('ECO dwarfs:', len(dwarfseleco))
###############################################################################

#S06 sample of RESOLVE and ECO from MPA-JHU
s06res = pd.read_csv("RESOLVE_full_bpt1snr5_dext_jhu.csv")
dwarfs06res = s06res.loc[s06res.logmstar < 9.5]
s06eco = pd.read_csv("ECO_full_bpt1snr5_dext_jhu.csv")
dwarfs06eco = s06eco.loc[s06eco.logmstar < 9.5]
print('Number of Galaxies in MPA-JHU S06 sample:')
print('RESOLVE :', len(s06res))
print('RESOLVE dwarfs:', len(dwarfs06res))
print('Number of Galaxies:')
print('ECO :', len(s06eco))
print('ECO dwarfs:', len(dwarfs06eco))

#S06 sample of RESOLVE and ECO from Portsmouth
s06res = pd.read_csv("RESOLVE_full_bpt1snr5_dext_port.csv")
dwarfs06res = s06res.loc[s06res.logmstar < 9.5]
s06eco = pd.read_csv("ECO_full_bpt1snr5_dext_port.csv")
dwarfs06eco = s06eco.loc[s06eco.logmstar < 9.5]
print('Number of Galaxies in Portsmouth S06 sample:')
print('RESOLVE :', len(s06res))
print('RESOLVE dwarfs:', len(dwarfs06res))
print('Number of Galaxies:')
print('ECO :', len(s06eco))
print('ECO dwarfs:', len(dwarfs06eco))

#S06 sample of RESOLVE and ECO from NSA
s06res = pd.read_csv("RESOLVE_full_bpt1snr5_dext_nsa.csv")
dwarfs06res = s06res.loc[s06res.logmstar < 9.5]
s06eco = pd.read_csv("ECO_full_bpt1snr5_dext_nsa.csv")
dwarfs06eco = s06eco.loc[s06eco.logmstar < 9.5]
print('Number of Galaxies in NSA S06 sample:')
print('RESOLVE :', len(s06res))
print('RESOLVE dwarfs:', len(dwarfs06res))
print('Number of Galaxies:')
print('ECO :', len(s06eco))
print('ECO dwarfs:', len(dwarfs06eco))
