# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 11:46:56 2021

@author: mugdhapolimera
"""


# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 09:41:34 2020

@author: mugdhapolimera

Select Gemini observing sample to quantify starburst-AGN-merger connection in
dwarfs

1. Compact 
2. Extended
3. Starbursting
4. Moderately SFing
5. Strong Emission lines?
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io.idl import readsav
import matplotlib as mpl
label_size = 12
mpl.rcParams['xtick.labelsize'] = label_size 
mpl.rcParams['ytick.labelsize'] = label_size 


threshold = 10
resfull = pd.read_csv("ECO_inobssample.csv")
resfull.index = resfull.name
#inspring = (resfull.radeg > 8.75*15.) & (resfull.radeg < 15.75*15.)
dwarf = (resfull.logmstar < threshold)

rescat = readsav("eco_wresa_032918.dat")
rescat['name'] = rescat['econames']
#rescatphot = readsav("resolvecatalogphot.dat")
resphot = pd.read_csv("ECO_barysample_photometrics.csv")
resphot.index = resphot.name
resndx = [x for x in range(len(rescat.name)) if rescat.name[x] in list(resfull.name)]
resphotndx = [x for x in range(len(resphot.name)) if resphot.name.iloc[x] in list(resfull.name)]

resfull['logmbary'] = np.log10(10**resfull.logmstar + 10**resfull.logmgas)
resfull['r50'] = rescat.r50[resndx]
resfull['mur90'] = rescat.mur90[resndx]
resfull['DEL_MU50'] = resphot.DEL_MU50[resphotndx]
resfull['DEL_SFR'] = resphot.DEL_SFR[resphotndx]
resfull['sfr_nuv_wise'] = resphot.SFR[resphotndx]
resfull['logssfr'] = np.log10(resfull['sfr_nuv_wise']) + 9 - resfull.logmstar
resfull['logsfr'] = np.log10(resfull['sfr_nuv_wise']) + 9 
#resfull['BN'] = resphot.BN[resphotndx].astype('bool')

#gem_red = ['rf0266','rf0284','rf0250','rf0045','rf0272','rf0370']
#sami = list(pd.read_csv("SAMI_RESOLVE.csv").name_1)
#manga = list(pd.read_csv("MANGA_RESOLVE.csv").name_1)
#trans= ['rs1036','rs0022','rs1287','rs0014','rs0320','rs0537']
#resfull = resfull[resfull.del_mu > 0]
#resfull.loc['rs1456']['BN'] = 0

res = resfull[dwarf & (resfull.logmbary > 9.2)] #inspring & inspring & 
dwarf = (res.logmstar < threshold)

resflag = pd.read_csv("eco_emlineclass_dext_snr5_jhu.csv")
resflag.index = resflag.galname
selagn = list(resflag.galname[resflag.sftoagn | resflag.agntosf | \
                              resflag.composite | resflag.defagn])
sels = list(resflag.galname)
ressel = resfull.loc[sels]#[inspring]
seldwarf = list(ressel[ressel.logmstar < 9.5].name)
seldwarfagn = list(np.intersect1d(selagn, seldwarf))
selim = list(ressel[(ressel.logmstar > 9.5) & (ressel.logmstar < 10)].name)
selimagn = list(np.intersect1d(selagn, selim))
#seldwarf = list(ressel[ressel.logmstar < threshold].name)
#sel_trans = np.intersect1d(sels,trans)

#barro = pd.read_csv("Barro_inobssample.csv")
#barro.index = barro.name 
#barro = barro.loc[sels]

#sel_sami = np.intersect1d(sels,sami)
#sel_manga = np.intersect1d(sels,manga)
#sel_gem = np.intersect1d(sels,gem_red)

highsfr = list(ressel.name[ressel.logssfr > -0.5])
bounds = [7.5,8.0,8.5,9., 9.5, 10.0]#np.arange(8.4,10.1,0.2)
#cmap = plt.get_cmap('Spectral',6)#int(np.max(r)/0.5))
cmap = mpl.colors.ListedColormap(['navy','navy','blue','mediumturquoise',\
                                  'darkgreen'])
boundaries = bounds#np.array([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5])#, 4.0])#, 4.5, 5.0])
norm = mpl.colors.BoundaryNorm(boundaries, cmap.N, clip=True)

#plt.figure()
#plt.ylabel('log($\Delta$SFR$_{MS}$)')
#plt.xlabel('$\Delta \mu_\Delta$')
#plt.plot(barro.DEL_MUDEL.loc[sels], \
#         barro.DEL_SFR.loc[sels], '.', color = 'gray', alpha = 0.5,
#         label = 'SEL')
#plt.plot(barro.DEL_MUDEL.loc[seldwarf], \
#         barro.DEL_SFR.loc[seldwarf], '.', color = 'blue', label = 'Dwarf SELs')
##plt.plot(res.DEL_MUDEL[res.BN], res.DEL_SFR[res.BN], 'o', color = 'blue', \
##         label = 'Blue Nuggets')
#plt.plot(barro.DEL_MUDEL.loc[seldwarfagn], \
#         barro.DEL_SFR.loc[seldwarfagn], 'o', color = 'blue', ms = 10, 
#         label = 'Optical Dwarf SEL AGN')
#plt.plot(barro.DEL_MUDEL.loc[selim], \
#         barro.DEL_SFR.loc[selim], '.', color = 'green', 
#         label = 'Intermediate Mass SEL')
#plt.plot(barro.DEL_MUDEL.loc[selimagn], \
#         barro.DEL_SFR.loc[selimagn], 'o', color = 'green', ms = 10,
#         label = 'Intermediate Mass SEL AGN')
##plt.plot(barro.DEL_MUDEL.loc[highsfr], \
##         barro.DEL_SFR.loc[highsfr], 'o', color = 'lime', 
##         label = 'High SSFR SEL')
#plt.legend(fontsize = 10)

#temp = rs1081,'rs0337',''rs0933'
#dwarfagntargets = ['rs0010','rs1047','rs1143','rs0909','rs1038','rs0472'] #'rs0421',
#dwarfsftargets = ['rs1163','rs1014','rs1200','rs0432',\
#                  'rs1091','rs1178']

#imagntargets = ['rs1150','rs1004','rs0978','rs1036']
#imsftargets = ['rs0087','rs0237','rs0096','rs0320']

anno = 0
#transition = ['rs0022', 'rs0320', 'rs1036']
#targets = dwarfagntargets + dwarfsftargets + imagntargets + imsftargets #+ list(sel_trans)
pottarg = np.intersect1d(resfull.name,seldwarf+selim)
fig,ax = plt.subplots()
#plt.plot(barro.DEL_MU50.loc[targets], \
#         barro.DEL_SFR.loc[targets], 'o', ms = 15, mec = 'r', mfc = 'none',
#         mew = 2, label = 'This Proposal Targets')
#plt.legend(fontsize = 15)
cax = ax.scatter(resfull.DEL_MU50.loc[seldwarf+selim],
                 resfull.DEL_SFR.loc[seldwarf+selim], marker = 'o', s=10, 
                 c = ressel.logmstar.loc[seldwarf+selim], cmap = cmap, norm = norm,
                 label = 'SEL')
cax2 = ax.scatter(resfull.DEL_MU50.loc[seldwarfagn+selimagn], 
                  resfull.DEL_SFR.loc[seldwarfagn+selimagn], marker = 'o', s=100, 
                 c = ressel.logmstar.loc[seldwarfagn+selimagn], cmap = cmap, norm= norm)
#cb2 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,norm=norm)
#cb2.set_label('stellar mass')
fig.colorbar(cax,extend = 'min', ticks = bounds,boundaries = bounds)
plt.ylabel('log($\Delta$SFR$_{MS}$)',fontsize = 15)
plt.xlabel('$\Delta \Sigma_e$',fontsize = 15)
#plt.plot(barro.DEL_MU50.loc[sels], \
#         barro.DEL_SFR.loc[sels], '.', color = 'gray', alpha = 0.5,
#         label = 'SEL')
#plt.plot(barro.DEL_MU50.loc[seldwarf], \
#         barro.DEL_SFR.loc[seldwarf], '.', color = 'blue', ms = 8, label = 'Dwarf SELs')
##plt.plot(res.DEL_MU50[res.BN], res.DEL_SFR[res.BN], 'o', color = 'blue', \
##         label = 'Blue Nuggets')
#plt.plot(barro.DEL_MU50.loc[seldwarfagn], \
#         barro.DEL_SFR.loc[seldwarfagn], 'o', color = 'blue', ms = 10, 
#         label = 'Optical Dwarf SEL AGN')
#plt.plot(barro.DEL_MU50.loc[selim], \
#         barro.DEL_SFR.loc[selim], '.', color = 'green', ms = 8, 
#         label = 'Intermediate Mass SEL')
#plt.plot(barro.DEL_MU50.loc[selimagn], \
#         barro.DEL_SFR.loc[selimagn], 'o', color = 'green', ms = 10,
#         label = 'Intermediate Mass SEL AGN')
##plt.plot(barro.DEL_MU50.loc[sel_manga], \
##         barro.DEL_SFR.loc[sel_manga], 'o', ms = 11, mec = 'gray', mfc = 'none',
##         mew = 2, label = 'SAMI/MaNGA/Gemini')
#plt.plot(barro.DEL_MU50.loc[targets], \
#         barro.DEL_SFR.loc[targets], 'o', ms = 11, mec = 'r', mfc = 'none',
#         mew = 2, label = 'This Proposal Targets')
#plt.plot(barro.DEL_MU50.loc[trans], \
#         barro.DEL_SFR.loc[trans], 'o', ms = 11, mec = 'magenta', mfc = 'none',
#         mew = 2, label = 'Transitional nuggets')
if anno ==1:
    for i, txt in enumerate(pottarg):
        plt.annotate(txt, (resfull.DEL_MU50.loc[pottarg[i]], \
         resfull.DEL_SFR.loc[pottarg[i]]), fontsize = 8)
#plt.plot(barro.DEL_MU50.loc[highsfr], \
#         barro.DEL_SFR.loc[highsfr], 'o', color = 'lime', 
#         label = 'High SSFR SEL')

#plt.plot(barro.DEL_MU50.loc[sel_sami], \
#         barro.DEL_SFR.loc[sel_sami], 'o', ms = 10, mec = 'gray', mfc = 'none',
#         mew = 2)
#plt.plot(barro.DEL_MU50.loc[sel_gem], \
#         barro.DEL_SFR.loc[sel_gem], 'o', ms = 10, mec = 'gray', mfc = 'none', 
#         mew = 2)


#plt.figure()
#plt.ylabel('log(M$_{bary}$/M$_\odot$)')
#plt.xlabel('log(M$_{halo}$/M$_\odot$)')
#plt.plot(ressel.logmh.loc[sels], \
#         ressel.logmbary.loc[sels], '.', color = 'gray', label = 'SEL')
#plt.plot(ressel.logmh.loc[seldwarf], \
#         ressel.logmbary.loc[seldwarf], '.', color = 'blue', label = 'Dwarf SELs')
##plt.plot(res.logmh[res.BN], res.logmbary[res.BN], 'o', color = 'blue', \
##         label = 'Blue Nuggets')
#plt.plot(ressel.logmh.loc[seldwarfagn], \
#         ressel.logmbary.loc[seldwarfagn], 'o',ms = 10, color = 'blue', 
#         label = 'Optical Dwarf SEL AGN')
#plt.plot(ressel.logmh.loc[selim], \
#         ressel.logmbary.loc[selim], '.', color = 'green', 
#         label = 'Intermediate Mass SEL')
#plt.plot(ressel.logmh.loc[selimagn], \
#         ressel.logmbary.loc[selimagn], 'o', ms = 10,color = 'green', 
#         label = 'Intermediate Mass SEL AGN')
#plt.plot(ressel.logmh.loc[targets], ressel.logmbary.loc[targets], \
#         'o', ms = 11, mec = 'r', mfc = 'none', mew = 2,
#         label = 'This Proposal Targets')
##plt.axhline(y=9.5, xmin=0.0, xmax=1.0, color='k')
##plt.axhline(y=10, xmin=0.0, xmax=1.0, color='k', ls = '--')
#plt.axvline(x=11.5, ymin=0.0, ymax=1.0, color='k')
#plt.axvline(x=12.0, ymin=0.0, ymax=1.0, color='k', ls = '--')
#plt.legend(fontsize = 10)

ressel['g_s'] = ressel.logmgas - ressel.logmstar
fig,ax = plt.subplots()
#plt.plot(ressel.logmstar.loc[targets], ressel.g_s.loc[targets], \
#         'o', ms = 11, mec = 'r', mfc = 'none', mew = 2,
#         label = 'This Proposal Targets')
plt.axhline(y=0, xmin=0.0, xmax=1.0, color='k', ls = '--')
#plt.axvline(x=9.5, ymin=0.0, ymax=1.0, color='k',ls = '--')
#plt.legend(fontsize = 10)
cax = ax.scatter(ressel.logmstar.loc[seldwarf+selim], ressel.g_s.loc[seldwarf+selim], marker = 'o', s=10, 
                 c = ressel.logmstar.loc[seldwarf+selim], cmap = cmap, norm = norm,
                 label = 'SEL')
cax2 = ax.scatter(ressel.logmstar.loc[seldwarfagn+selimagn], ressel.g_s.loc[seldwarfagn+selimagn], marker = 'o', s=100, 
                 c = ressel.logmstar.loc[seldwarfagn+selimagn], cmap = cmap, norm= norm)

plt.ylabel('log(M$_{gas}$/M$_*$)')
plt.xlabel('log(M$_{*}$/M$_\odot$)')
if anno == 1:
    for i, txt in enumerate(pottarg):
        plt.annotate(txt, (ressel.logmstar.loc[pottarg[i]], \
             ressel.g_s.loc[pottarg[i]]), fontsize = 8)

#plt.plot(ressel.logmstar.loc[sels], \
#         ressel.g_s.loc[sels], '.', color = 'gray', label = 'SEL')
#plt.plot(ressel.logmstar.loc[seldwarf], \
#         ressel.g_s.loc[seldwarf], '.', color = 'blue', label = 'Dwarf SELs')
##plt.plot(res.logmstar[res.BN], res.g_s[res.BN], 'o', color = 'blue', \
##         label = 'Blue Nuggets')
#plt.plot(ressel.logmstar.loc[seldwarfagn], \
#         ressel.g_s.loc[seldwarfagn], 'o',ms = 10, color = 'blue', 
#         label = 'Optical Dwarf SEL AGN')
#plt.plot(ressel.logmstar.loc[selim], \
#         ressel.g_s.loc[selim], '.', color = 'green', 
#         label = 'Intermediate Mass SEL')
#plt.plot(ressel.logmstar.loc[selimagn], \
#         ressel.g_s.loc[selimagn], 'o', ms = 10,color = 'green', 
#         label = 'Intermediate Mass SEL AGN')

#bins = np.arange(7.7,11.5,0.2)
#plt.figure()
#plt.hist(ressel.logmstar.loc[sels], bins = bins, color = 'gray',label = 'SEL', density = True)
#plt.hist(ressel.logmstar.loc[seldwarf], color = 'blue', label = 'Dwarf SELs', 
#         density = True, histtype = 'step', linewidth = 2)
#plt.hist(ressel.logmstar.loc[seldwarfagn], bins = bins, color = 'darkblue', 
#         label = 'Optical Dwarf SEL AGN', density = True, histtype = 'step', 
#         hatch = '\\', linewidth = 2)
#plt.hist(ressel.logmstar.loc[selim],bins = bins, color = 'limegreen',
#         histtype = 'step', linewidth = 2,
#         label = 'Intermediate Mass SEL', density = True)
#plt.hist(ressel.logmstar.loc[selimagn], bins = bins, color = 'darkgreen', 
#         histtype = 'step', hatch = '/', linewidth = 2,
#         label = 'Intermediate Mass SEL AGN', density = True)
#plt.hist(ressel.logmstar.loc[targets], bins = bins, color = 'red', 
#         histtype = 'step', linewidth = 2,
#         label = 'This Proposal Targets', density = True)
#plt.xlabel('log(M$_{*}$/M$_\odot$)')
#plt.ylabel('Relative frequency')
#plt.legend()
#
#bins = np.arange(10.5,14.3,0.2)
#plt.figure()
#plt.hist(ressel.logmh.loc[sels], bins = bins, color = 'gray',label = 'SEL', density = True)
#plt.hist(ressel.logmh.loc[seldwarf], color = 'blue', label = 'Dwarf SELs', 
#         density = True, histtype = 'step', linewidth = 2)
#plt.hist(ressel.logmh.loc[seldwarfagn], bins = bins, color = 'darkblue', 
#         label = 'Optical Dwarf SEL AGN', density = True, histtype = 'step', 
#         hatch = '\\', linewidth = 2)
#plt.hist(ressel.logmh.loc[selim],bins = bins, color = 'limegreen',
#         histtype = 'step', linewidth = 2,
#         label = 'Intermediate Mass SEL', density = True)
#plt.hist(ressel.logmh.loc[selimagn], bins = bins, color = 'darkgreen', 
#         histtype = 'step', hatch = '/', linewidth = 2,
#         label = 'Intermediate Mass SEL AGN', density = True)
#plt.hist(ressel.logmh.loc[targets], bins = bins, color = 'red', 
#         histtype = 'step', linewidth = 2,
#         label = 'This Proposal Targets', density = True)
#plt.xlabel('log(M$_{halo}$/M$_\odot$)')
#plt.ylabel('Relative Frequency')
#plt.legend()

ssfr_st = 10**(np.log10(ressel.sfr_nuv_wise) - ressel.logmstar)
ssfr_lt = 1/(1+(1/(ressel.meanfsmgr)))
fsmgr_st = 100*(10**6)*(ssfr_st)/(1-ssfr_st)/(0.1*1e9)
fsmgr_lt = ressel.meanfsmgr


fig,ax = plt.subplots()
#plt.plot(fsmgr_st.loc[targets], fsmgr_lt.loc[targets], \
#         'o', ms = 15, mec = 'r', mfc = 'none', mew = 2,
#         label = 'This Proposal Targets')
#plt.legend(fontsize = 15, loc = 'lower right')
cax = ax.scatter(fsmgr_st.loc[seldwarf+selim], fsmgr_lt.loc[seldwarf+selim], marker = 'o', s=10, 
                 c = ressel.logmstar.loc[seldwarf+selim], cmap = cmap, norm = norm,
                 label = 'SEL')
cax2 = ax.scatter(fsmgr_st.loc[seldwarfagn+selimagn], fsmgr_lt.loc[seldwarfagn+selimagn], marker = 'o', s=100, 
                 c = ressel.logmstar.loc[seldwarfagn+selimagn], cmap = cmap, norm= norm,
                 label = 'SEL')
#cb2 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,norm=norm)
#cb2.set_label('stellar mass')
fig.colorbar(cax,ax=[ax],extend = 'min', ticks = bounds,boundaries = bounds)
xaxis = np.arange(np.min(fsmgr_st)-0.05, np.max(fsmgr_st)+0.05,0.01)
yaxis = np.ones(len(xaxis))
#plt.plot(xaxis, yaxis, 'k--', lw = 3)
plt.xlim(np.min(fsmgr_st), np.max(fsmgr_st))
plt.ylim(np.min(fsmgr_lt), np.max(fsmgr_lt))
#plt.text(0.0005, 1.25, r'Stellar Mass Doubled in last Gyr', 
#             fontsize=14, color='k')
plt.yscale('log')
plt.xscale('log')
plt.xlabel(r'Short Term SFH', fontsize = 18)
# $\left(\frac{M_*(<100 Myr)}{M_*(>100 Myr)}\right)$',
plt.ylabel(r'Long Term SFH ',fontsize = 18)
fit = np.poly1d(np.polyfit(np.log10(fsmgr_st), \
                     np.log10(fsmgr_lt), 1))
y = np.log10(xaxis)*fit[1] + fit[0]
#plt.plot(xaxis,10**y,ls='--',c='k')
#$\left(\frac{M_*(<1 Gyr)}{M_*(>1 Gyr)}\right)$'
################################
#
# Define demarcation functions
#

#def o3hbcomposite(log_NII_HA):
#    return (0.61 / (log_NII_HA - 0.05)) + 1.3
#
#def o3hbmain(log_NII_HA):
#    return (0.61 / (log_NII_HA - 0.47)) + 1.19
#
#def o1hamain(log_OI_HA): #main line for OI/H-alpha from equation 3, Kewley 2006
#    return 1.33 + (0.73 / (log_OI_HA + 0.59))
#def o1hacrit(log_OI_HA): #boundary for OI/H-alpha
#    return -0.59
#def s2hamain(log_SII_HA): #main line for SII/H-alpha from equation 2, Kewley 2006
#    return 1.30 + (0.72 / (log_SII_HA - 0.32))

def n2hacompmin(log_NII_HA): #composite minimum line from equation 1, Kewley 2006
    return 1.3 + (0.61 / (log_NII_HA - 0.05))
def n2hamain(log_NII_HA): #main line for NII/H-alpha from equation 5, Kewley 2006
    return 1.19 + (0.61 / (log_NII_HA - 0.47))
def s2hamain(log_SII_HA): #main line for SII/H-alpha from equation 2, Kewley 2006
    return 1.30 + (0.72 / (log_SII_HA - 0.32))
def s2halinseyf(log_SII_HA): #liner/seyfert divider for SII/H-alpha
    return 0.76 + 1.89*log_SII_HA
def o1hamain(log_OI_HA): #main line for OI/H-alpha from equation 3, Kewley 2006
    return 1.33 + (0.73 / (log_OI_HA + 0.59))
def o1halinseyf(log_OI_HA): #liner/seyfert divider for OI/H-alpha
    return 1.3 + 1.18*log_OI_HA
def o1hacrit(log_OI_HA): #boundary for OI/H-alpha
    return -0.59

refn2ha = np.linspace(-3.0, 0.35)
refoiha = np.linspace(-2.5, -0.4)
refsiiha = np.linspace(-2, 0.3,100)
main_sii = s2hamain(refsiiha)
main_oi = o1hamain(refoiha)

f, (sp1, sp2, sp3) = plt.subplots(1,3, sharey = True)
#bounds = [7.8,8.4, 8.7, 9. , 9.3, 9.6, 10.0]#np.arange(8.4,10.1,0.2)
##cmap = plt.get_cmap('Spectral',6)#int(np.max(r)/0.5))
#cmap = mpl.colors.ListedColormap(['navy','navy','blue','lightsteelblue','mediumturquoise',\
#                                  'darkgreen'])
#boundaries = bounds#np.array([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5])#, 4.0])#, 4.5, 5.0])
#norm = mpl.colors.BoundaryNorm(boundaries, cmap.N, clip=True)
grp_color = ['blue','blue','green','green','red']
#grps = [seldwarf,seldwarfagn,selim,selimagn,targets]
grps = [seldwarf+selim,seldwarfagn+selimagn,targets]
grp_size = [10,100,15]
grp_marker = ['o','o','o']#,'o','o']
grp_label = ['Dwarf SELs','Dwarf SEL AGN','Proposal Targets']
a = 1
niiha = np.log10(ressel['nii_6584_flux']/ressel['h_alpha_flux'])
SII = ressel['sii_6731_flux'] + ressel['sii_6717_flux']
siiha = np.log10(SII/ressel['h_alpha_flux'])
oiha = np.log10(ressel['oi_6300_flux']/ressel['h_alpha_flux'])
oiiihb = np.log10(ressel['oiii_5007_flux']/ressel['h_beta_flux'])

for i in range(len(grps)):
    if i == 2:
        sp1.plot(niiha.loc[grps[i]],oiiihb.loc[grps[i]], 'o',c=grp_color[i], ms= grp_size[i],
                 mfc = 'none',mec = 'r',mew = 2, label = grp_label[i])
        sp2.plot(siiha.loc[grps[i]],oiiihb.loc[grps[i]], 'o',c=grp_color[i], ms= grp_size[i],mfc = 'none',mec = 'r',mew = 2, 
                    label = grp_label[i])
        sp3.plot(oiha.loc[grps[i]],oiiihb.loc[grps[i]], 'o',c=grp_color[i], ms= grp_size[i],mfc = 'none',mec = 'r',mew = 2)
    else: 
        sp1.scatter(niiha.loc[grps[i]],oiiihb.loc[grps[i]], marker = grp_marker[i], s= grp_size[i], 
                    alpha = a,c = ressel.logmstar.loc[grps[i]],
                    cmap = cmap, norm= norm, label = grp_label[i])
        sp2.scatter(siiha.loc[grps[i]],oiiihb.loc[grps[i]], marker = grp_marker[i], s= grp_size[i], 
                    alpha = a,c = ressel.logmstar.loc[grps[i]],
                    cmap = cmap, norm= norm)
        sp3.scatter(oiha.loc[grps[i]].loc[grps[i]],oiiihb.loc[grps[i]], marker = grp_marker[i], s= grp_size[i], 
                    alpha = a,c = ressel.logmstar.loc[grps[i]],
                    cmap = cmap, norm= norm)
#for i in range(len(grps)):
#    niiha = np.log10(ressel.loc[grps[i]]['nii_6584_flux']/ressel.loc[grps[i]]['h_alpha_flux'])
#    SII = ressel.loc[grps[i]]['sii_6731_flux'] + ressel.loc[grps[i]]['sii_6717_flux']
#    siiha = np.log10(SII/ressel.loc[grps[i]]['h_alpha_flux'])
#    oiha = np.log10(ressel.loc[grps[i]]['oi_6300_flux']/ressel.loc[grps[i]]['h_alpha_flux'])
#    oiiihb = np.log10(ressel.loc[grps[i]]['oiii_5007_flux']/ressel.loc[grps[i]]['h_beta_flux'])
#    
#    if i == 4:
#        sp1.plot(niiha,oiiihb, 'o',c=grp_color[i], ms= grp_size[i],mfc = 'none',mec = 'r',mew = 2)
#        sp2.plot(siiha,oiiihb, 'o',c=grp_color[i], ms= grp_size[i],mfc = 'none',mec = 'r',mew = 2, 
#                    label = grp_label[i])
#        sp3.plot(oiha,oiiihb, 'o',c=grp_color[i], ms= grp_size[i],mfc = 'none',mec = 'r',mew = 2)
#    else: 
#        sp1.scatter(niiha,oiiihb, c=grp_color[i],marker = grp_marker[i], s= grp_size[i], alpha = 0.5)
#        sp2.scatter(siiha,oiiihb, c=grp_color[i],marker = grp_marker[i], s= grp_size[i],  alpha = 0.5,
#                    label = grp_label[i])
#        sp3.scatter(oiha,oiiihb, c=grp_color[i],marker = grp_marker[i], s= grp_size[i], alpha = 0.5)
sp1.plot(refn2ha, n2hamain(refn2ha),'k',zorder = 0)
sp1.plot(refn2ha[refn2ha < 0], n2hacompmin(refn2ha[refn2ha < 0]),'k--',zorder = 0)
sp2.plot(refsiiha, s2hamain(refsiiha),'k',zorder = 0)
sp2.plot(refsiiha[refsiiha > -0.31], s2halinseyf(refsiiha[refsiiha > -0.31]),
              'k-.',zorder = 0)
    
sp2.set_xlabel(r'$\rm log([$SII$]$ 6717 + 6731/H$\alpha)$',fontsize=22)
sp3.plot(refoiha[refoiha < -0.7], o1hamain(refoiha[refoiha < -0.7]),'k', 
         zorder = 0)
sp3.plot(refoiha[refoiha > -1.13], o1halinseyf(refoiha[refoiha > -1.13]),
                           'k-.',zorder = 0)

sp3.set_xlabel(r'$\rm log ([$OI$]$ 6300/H$\alpha$)',fontsize=22)
sp1.set_xlabel(r'$\rm log([$N II$]$ 6584 / H$\alpha)$',fontsize=22)
sp1.set_ylabel(r'$\rm log([$O III$]$ 5007 / H$\beta$)',fontsize=22)
#sp1.legend(loc = 'lower left', fontsize = 12)
xmin, xmax = -2.00001, 0.30001
ymin, ymax = -1.5, 1.5
sp1.set_xlim(xmin,xmax)
sp1.set_ylim(ymin,ymax)
xmin, xmax = -2.00001, 0.50001
sp2.set_xlim(xmin,xmax)
xmin, xmax = -2.50001, -0.40001
sp3.set_xlim(xmin,xmax)
if anno ==1:
    for i, txt in enumerate(pottarg):
        sp1.annotate(txt, (niiha.loc[pottarg[i]], \
         oiiihb.loc[pottarg[i]]), fontsize = 8)
#selndx = [x for x in range(len(rescat.name)) if rescat.name[x] in sels]
#targndx = [x for x in range(len(rescat.name)) if rescat.name[x] in targets]
#ressel['rmag'] = 0
#ressel.rmag.loc[sels] = np.array(resphotcat.rmag[selndx], dtype = np.float)        
#targ_table = pd.DataFrame({'Name': np.array(targets),
#                           'RA': ressel.radeg.loc[targets],
#                           'Dec': ressel.dedeg.loc[targets],
#                           'r': ressel.rmag.loc[targets],
#                           'r_sys': ['AB']*len(targets)})
#targ_table.to_csv("Gemini2020A_targets.csv",index = False)