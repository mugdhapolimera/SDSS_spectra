# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 08:18:27 2019

@author: mugdhapolimera
"""

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_filter_new.pkl'
df = pd.read_pickle(inputfile)

ax = plt.subplot(121)
ax.plot(df.Flux_Ha_6562, df.h_alpha_flux,'o')
ax.plot(np.linspace(0,30000),np.linspace(0,30000), 'r')
ax.set_xlim(0,20000)
ax.set_ylim(0,20000)
ax.set_xlabel(r'H$\alpha$ Portsmouth',size = 14)
ax.set_ylabel(r'H$\alpha$ MPA-JHU',size = 14)
ax = plt.subplot(122)
ax.set_xlim(0,20000)
ax.set_ylim(-20000,20000)
ax.plot(df.Flux_Ha_6562,df.Flux_Ha_6562-df.h_alpha_flux,'o')
ax.set_xlabel(r'H$\alpha$ Portsmouth',size = 14)
ax.set_ylabel(r'H$\alpha$ Residuals (Porstmouth - MPAJHU)',size = 14)
ax.plot(np.linspace(0,30000),0*np.linspace(0,30000), 'r')

ax = plt.figure()
ax = plt.subplot(121)
ax.plot(df.Flux_Hb_4861, df.h_beta_flux,'o')
ax.plot(np.linspace(0,10500),np.linspace(0,10500), 'r')
ax.set_xlim(0,6000)
ax.set_ylim(0,6000)
ax.set_xlabel(r'H$\beta$ Portsmouth',size = 14)
ax.set_ylabel(r'H$\beta$ MPA-JHU',size = 14)
ax = plt.subplot(122)
ax.set_xlim(0,6000)
ax.set_ylim(-6000,6000)
ax.plot(df.Flux_Hb_4861,df.Flux_Hb_4861-df.h_beta_flux,'o')
ax.set_xlabel(r'H$\beta$ Portsmouth',size = 14)
ax.set_ylabel(r'H$\beta$ Residuals (Porstmouth - MPAJHU)',size = 14)
ax.plot(np.linspace(0,6000),0*np.linspace(0,6000), 'r')

ax = plt.figure()
ax = plt.subplot(121)
ax.plot(df.Flux_Ha_6562/df.Flux_Hb_4861, df.h_alpha_flux/df.h_beta_flux,'o')
ax.plot(np.linspace(0,45),np.linspace(0,45), 'r')
ax.set_xlabel(r'H$\alpha$/H$\beta$ Portsmouth',size = 14)
ax.set_ylabel(r'H$\alpha$/H$\beta$ MPA-JHU',size = 14)
ax.set_xlim(0,40)
ax.set_ylim(0,40)
ax = plt.subplot(122)
ax.set_xlim(0,40)
ax.set_ylim(-40,40)
ax.plot(df.Flux_Ha_6562/df.Flux_Hb_4861,
        (df.Flux_Ha_6562/df.Flux_Hb_4861)-(df.h_alpha_flux/df.h_beta_flux),'o')
ax.set_xlabel(r'H$\alpha$/H$\beta$ Portsmouth',size = 14)
ax.set_ylabel(r'H$\alpha$/H$\beta$ Residuals (Porstmouth - MPAJHU)',size = 14)
ax.plot(np.linspace(0,45),0*np.linspace(0,45), 'r')
