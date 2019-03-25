# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 20:56:19 2018

@author: mugdhapolimera
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 20:56:19 2018

@author: mugdhapolimera
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#rdat = pd.read_pickle("C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_bpt1filter.pkl")
rdat = pd.read_pickle("C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_bpt1_filter.pkl")
# a simple Te method calib. is the NII/Halpha PP04, for ???2.5 < N2 < ???0.3
# Pettini & Pagel 2004

N2 = np.log10(rdat.nii_6584_flux/rdat.h_alpha_flux)

#selgd<-which(N2>(-2.5) & N2<(-0.3))
selbd = ((N2<-2.5) | (N2>-0.3))
PP04_12logOH = 9.37+2.03*N2+1.26*N2**2+0.32*N2**3
rep = np.isnan(PP04_12logOH)
PP04_12logOH[rep] = (-99.)
PP04_12logOH[selbd] = (-99.)

# alternative PP04 calib using O3, hbeta, and N, Halpha
O3N2 = np.log10((rdat.oiii_5007_flux/rdat.h_beta_flux)/(rdat.nii_6584_flux/rdat.h_alpha_flux))

PP04_12logOH_v2 = 8.73-0.32*O3N2

rep = np.isnan(PP04_12logOH_v2)
PP04_12logOH_v2[rep] = (-99.)
selbd = (O3N2 > 1.9) # supposedly this calibration is bad in this regime, don't use those values
PP04_12logOH_v2[selbd] = (-99.)

#plt.plot(PP04_12logOH, PP04_12logOH_v2,'o')
#plt.xlabel("[NII]/Halpha Z")
#plt.ylabel("[OIII]/Hbeta &[NII]/Halpha Z")
#plt.xlim(7.5,10.5)
#plt.ylim(7.5,10.5)
#plt.plot(np.arange(7.5,11),np.arange(7.5,11), "r")
#


#out<-cbind(as.vector(rdat[,"NAME"]), PP04_12logOH, PP04_12logOH_v2)

#write.csv(out, "RES_PP04Zest.csv", row.names=FALSE, quote=FALSE)



#iziout = np.genfromtxt("C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_bpass_bpt1filter_SEL_newgrid.txt", dtype = None, names = ["name", 
#                                                                                         "Estimate", "err_up", "err_down"])
#iziout = np.genfromtxt("C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_bpass_bpt1filter_SEL_grovesdep.txt", dtype = None, names = ["name", 
#                                                                                         "Estimate", "err_up", "err_down"])
iziout = np.genfromtxt("C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_nicholls_bpt1filter_SEL.txt", dtype = None, names = ["name", 
                                                                                         "Estimate", "err_up", "err_down"])

iziout2 = np.genfromtxt("C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_bpass_bpt1filter_SEL_groves_jen.txt", dtype = None, names = ["name", 
                                                                                         "Estimate", "err_down", "err_up"])
matchi = np.arange(913)
#np.where(i for i in range(len(iziout)) if rdat.name[i] in iziout["name"])
#PP04_12logOH

#png(file="ResZcomp_PP04.png", width=6, height=6, units="in", res=200)
#gdpts= ~(np.isnan(matchi))
x_solar = np.linspace(7.4, 9)
y_solar_40 = 8.76+np.log10(0.4)*np.ones(len(x_solar))
y_solar_30 = 8.76+np.log10(0.3)*np.ones(len(x_solar))
y_solar_20 = 8.76+np.log10(0.2)*np.ones(len(x_solar))

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax2 = ax1.twiny()
ax1.plot(PP04_12logOH, iziout["Estimate"], 'ko', alpha = 0.25)
ax1.set_xlim(7.4,9)
ax1.set_ylim(7.4,9)
ax1.plot(x_solar, y_solar_40,'b')
ax1.plot(x_solar, y_solar_30,'b')
ax1.plot(x_solar, y_solar_20,'b')
ax1.plot(y_solar_40,x_solar,'b')
ax1.plot(y_solar_30,x_solar,'b')
ax1.plot(y_solar_20,x_solar,'b')

ax1.set_title("NB(BPASS+Nicholls) vs PP04", y = 1.1)
ax1.set_ylabel("12 + log(O/H)  [NB + L10]")
ax1.set_xlabel("12 + log(O/H)  [PP04]")
ax2.set_xlabel("N2 (NII/H-alpha)")
ax1.plot(np.arange(7.4,11),np.arange(7.4,11), "r")
#ax1.plot(np.arange(-1.3,0.5),np.arange(-1.3,0.5), "r")
ax2.set_xticks([0.0,0.2,0.4,0.6,0.8,1.0, 1.2])
float_formatter = lambda x: "%.2f" % x
xticks = np.array([7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0])
N2 = 1.754*xticks - 15.614
N2_label = ["%.2f" % z for z in N2]
<<<<<<< HEAD
ax2.set_xticklabels(N2_label)

fig = plt.figure()
ax1 = fig.add_subplot(111)
=======
N2_20 = 1.754*(np.log10(0.2)+8.76) - 15.614
N2_30 = 1.754*(np.log10(0.3)+8.76) - 15.614
N2_40 = 1.754*(np.log10(0.4)+8.76) - 15.614

fig = plt.figure()
ax1 = fig.add_subplot(122)
ax2 = ax1.twiny()
ax2.set_xticklabels(N2_label)
ax1.plot(np.arange(7.8,11),np.arange(7.8,11), "r", linewidth = 3)
ax1.plot(PP04_12logOH, iziout["Estimate"], 'bo', alpha = 0.25)
ax1.set_xlim(7.8,9)
ax1.set_ylim(7.8,9)
ax1.plot(x_solar, y_solar_40,'k-.', linewidth = 2)
ax1.plot(x_solar, y_solar_30,'k-.', linewidth = 2)
ax1.plot(x_solar, y_solar_20,'k-.', linewidth = 2)
ax1.text(7.825, y_solar_20[0]+0.025, r'0.2 $Z_{\odot}$', fontsize=14, color='k')
ax1.text(7.825, y_solar_30[0]+0.025, r'0.3 $Z_{\odot}$', fontsize=14, color='k')
ax1.text(7.825, y_solar_40[0]+0.025, r'0.4 $Z_{\odot}$', fontsize=14, color='k')
ax1.plot(y_solar_40,x_solar,'k-.', linewidth = 2)
ax1.plot(y_solar_30,x_solar,'k-.', linewidth = 2)
ax1.plot(y_solar_20,x_solar,'k-.', linewidth = 2)
ax1.text(y_solar_20[0]-0.05, 8.9, r'[NII]/H$\alpha$ = ' + "%0.2f" % N2_20,
         fontsize=14, color='k', rotation = 'vertical')
ax1.text(y_solar_30[0]-0.05, 8.9, r'[NII]/H$\alpha$ = ' + "%0.2f" % N2_30,
         fontsize=14, color='k', rotation = 'vertical')
ax1.text(y_solar_40[0]-0.05, 8.9, r'[NII]/H$\alpha$ = ' + "%0.2f" % N2_40,
         fontsize=14, color='k', rotation = 'vertical')
#ax1.set_title("NB(BPASS+Nicholls) vs PP04", y = 1.1)
ax1.set_ylabel("12 + log(O/H)  (using Cloudy/BPASS modelling)", size = 15)
ax1.set_xlabel(r"12 + log(O/H)  (using [NII]/H$\alpha$)", size = 15)
ax2.set_xlabel(r"[NII]/H$\alpha$", size = 15)
#ax1.plot(np.arange(-1.3,0.5),np.arange(-1.3,0.5), "r")
ax2.set_xticks([0.0,0.2,0.4,0.6,0.8,1.0,1.2])
ax2.set_xticklabels(N2_label)
size = fig.get_size_inches()*fig.dpi 
ax1.set_title('(b)', y = -0.15)

ax1 = fig.add_subplot(121)
>>>>>>> e822faa7233d07eaed16d018b6ba44eed4db7eaf
ax2 = ax1.twiny()
ax1.plot(PP04_12logOH, iziout2["Estimate"], 'bo', alpha = 0.25)
<<<<<<< HEAD
ax1.set_xlim(7.4,9)
ax1.set_ylim(7.4,9)
ax1.plot(x_solar, y_solar_40,'k')
ax1.plot(x_solar, y_solar_30,'k')
ax1.plot(x_solar, y_solar_20,'k')
ax1.plot(y_solar_40,x_solar,'k')
ax1.plot(y_solar_30,x_solar,'k')
ax1.plot(y_solar_20,x_solar,'k')
ax1.set_title("NB(BPASS+Groves+No N dep) vs PP04", y=1.1)
ax1.set_ylabel("12 + log(O/H)  [NB + Chris-STB99]")
ax1.set_xlabel("12 + log(O/H)  [PP04]")
ax2.set_xlabel("N2 (NII/H-alpha)")
ax1.plot(np.arange(7.4,11),np.arange(7.4,11), "r")
=======
ax1.set_xlim(7.8,9)
ax1.set_ylim(7.8,9)
ax1.plot(x_solar, y_solar_40,'k-.', linewidth = 2)
ax1.plot(x_solar, y_solar_30,'k-.', linewidth = 2)
ax1.plot(x_solar, y_solar_20,'k-.', linewidth = 2)
ax1.text(7.825, y_solar_20[0]+0.025, r'0.2 $Z_{\odot}$', fontsize=14, color='k')
ax1.text(7.825, y_solar_30[0]+0.025, r'0.3 $Z_{\odot}$', fontsize=14, color='k')
ax1.text(7.825, y_solar_40[0]+0.025, r'0.4 $Z_{\odot}$', fontsize=14, color='k')
ax1.plot(y_solar_40,x_solar,'k-.', linewidth = 2)
ax1.plot(y_solar_30,x_solar,'k-.', linewidth = 2)
ax1.plot(y_solar_20,x_solar,'k-.', linewidth = 2)
ax1.text(y_solar_20[0]-0.05, 8.9, r'[NII]/H$\alpha$ = ' + "%0.2f" % N2_20,
         fontsize=14, color='k', rotation = 'vertical')
ax1.text(y_solar_30[0]-0.05, 8.9, r'[NII]/H$\alpha$ = ' + "%0.2f" % N2_30,
         fontsize=14, color='k', rotation = 'vertical')
ax1.text(y_solar_40[0]-0.05, 8.9, r'[NII]/H$\alpha$ = ' + "%0.2f" % N2_40,
         fontsize=14, color='k', rotation = 'vertical')
#ax1.set_title("NB(Levesque10) vs PP04", y=1.1)
#ax1.set_ylabel("12 + log(O/H)  (using Levesque 2010)", size = 15)
ax1.set_ylabel("12 + log(O/H)  (using grid with Hydrogen Density)", size = 15)
ax1.set_xlabel(r"12 + log(O/H)  (using [NII]/H$\alpha$)", size = 15)
ax2.set_xlabel(r"[NII]/H$\alpha$", size = 15)
>>>>>>> e822faa7233d07eaed16d018b6ba44eed4db7eaf
#ax1.plot(np.arange(7.7,11),np.arange(-1.3,0.5), "r")
ax2.set_xticklabels(N2_label)
