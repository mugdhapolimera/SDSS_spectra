# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 12:08:34 2022

@author: mugdhapolimera
"""

import numpy as np
from scipy.io.idl import readsav
import pandas as pd
import matplotlib.pyplot as plt

res = pd.read_csv('RESOLVE_WISE_good.csv')
res.index= res.name
full = pd.read_csv('ECO+RESOLVE_WISE_good.csv')
full.index = full.name
rescommon = np.intersect1d(res.name, full.name)
print(rescommon)

ecosel = pd.read_csv('ECO_SEL_WISE_good.csv')
ecosel.index= ecosel.name
ecoselcommon = np.intersect1d(ecosel['name'], full.name)
print(ecoselcommon)

mw1res = ecosel.loc[ecoselcommon]['mw1'] - full.loc[ecoselcommon]['mw1']
mw1res_err = np.sqrt(ecosel.loc[ecoselcommon]['emw1']**2 + full.loc[ecoselcommon]['emw1']**2)


plt.figure()
plt.plot(ecosel.loc[ecoselcommon]['mw1'], ecosel.loc[ecoselcommon]['mw1'] - \
         full.loc[ecoselcommon]['mw1'], 'o')
plt.errorbar(ecosel.loc[ecoselcommon]['mw1'], mw1res, yerr = mw1res_err, fmt = 'o')

diff = mw1res.index[np.abs(mw1res) > 0.01]

diffjushr = [ecodat['names'][i] for i in range(len(ecodat['names'])) \
             if ecodat['econames'][i] in diff]


resecoselcommon = np.intersect1d(res.econame, ecosel.name)
reseco = res[res.econame != 'notineco']
reseco.index = reseco.econame
mw1res = reseco.loc[resecoselcommon]['mw1'] - ecosel.loc[resecoselcommon]['mw1']
mw1res_err = np.sqrt(reseco.loc[resecoselcommon]['emw1']**2 + \
                     ecosel.loc[resecoselcommon]['emw1']**2)


plt.figure()
plt.errorbar(reseco.loc[ecoselcommon]['mw1'], mw1res, yerr = mw1res_err, fmt = 'o')

diff = mw1res.index[np.abs(mw1res) > 0.01]

#diffjushr = [ecodat['names'][i] for i in range(len(ecodat['names'])) \
#             if ecodat['econames'][i] in diff]


