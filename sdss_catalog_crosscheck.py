# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 11:11:34 2020

@author: mugdhapolimeradeg

SDSS catalog cross-check
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore") 
from scipy import stats
eco = 0
resolve = 1
full = 0
if eco: 
    eco = pd.read_csv('ECO_live22Oct2018.csv')
    eco.index = eco.name
    jhuflag = pd.read_csv('ECO/SEL/eco_emlineclass_dext_snr5_jhu.csv')
    jhuflag.index = jhuflag.galname
    jhuflagagn = jhuflag.sftoagn | jhuflag.agntosf | jhuflag.defagn | jhuflag.composite
    portflag = pd.read_csv('ECO/SEL/eco_emlineclass_dext_snr5_port.csv')
    portflag.index = portflag.galname
    portflagagn = portflag.sftoagn | portflag.agntosf | portflag.defagn | portflag.composite
    nsaflag = pd.read_csv('ECO/SEL/eco_emlineclass_dext_snr5_nsa.csv')
    nsaflag.index = nsaflag.galname
    nsaflagagn = nsaflag.sftoagn | nsaflag.agntosf | nsaflag.defagn | nsaflag.composite
    
    port = pd.read_csv('ECO/SEL/ECO_full_snr5_dext_port.csv')#[portflag.sftoagn]
    port.index = port.name
    jhu = pd.read_csv('ECO/SEL/ECO_full_snr5_dext_jhu.csv')#[jhuflag.sftoagn]
    jhu.index = jhu.name
    nsa = pd.read_csv('ECO/SEL/ECO_full_snr5_dext_nsa.csv')#.loc[nsaflag.index.values[nsaflag.sftoagn]]
    nsa.index = nsa.name
    master = pd.read_csv('ECO_snr5_master_bary.csv')
    master.index = master.name

if resolve: 
    eco = pd.read_csv('RESOLVE_liveMay2020.csv')
    eco = pd.read_csv('ECO_live22Oct2018.csv')
    eco.index = eco.name
    jhuflag = pd.read_csv('resolve_emlineclass_dext_snr5_jhu.csv')
    jhuflag.index = jhuflag.galname
    jhuflagagn = jhuflag.sftoagn | jhuflag.agntosf | jhuflag.defagn | jhuflag.composite
    
    portflag = pd.read_csv('resolve_emlineclass_dext_snr5_port.csv')
    portflag.index = portflag.galname
    portflagagn = portflag.sftoagn | portflag.agntosf | portflag.defagn | portflag.composite
    
    nsaflag = pd.read_csv('resolve_emlineclass_dext_snr5_nsa.csv')
    nsaflag.index = nsaflag.galname
    nsaflagagn = nsaflag.sftoagn | nsaflag.agntosf | nsaflag.defagn | nsaflag.composite
    jhuflag = pd.read_csv('ECO/SEL/eco_emlineclass_dext_snr5_jhu.csv')
    jhuflag.index = jhuflag.galname
    jhuflagagn = jhuflag.sftoagn | jhuflag.agntosf | jhuflag.defagn | jhuflag.composite
    portflag = pd.read_csv('ECO/SEL/eco_emlineclass_dext_snr5_port.csv')
    portflag.index = portflag.galname
    portflagagn = portflag.sftoagn | portflag.agntosf | portflag.defagn | portflag.composite
    nsaflag = pd.read_csv('ECO/SEL/eco_emlineclass_dext_snr5_nsa.csv')
    nsaflag.index = nsaflag.galname
    nsaflagagn = nsaflag.sftoagn | nsaflag.agntosf | nsaflag.defagn | nsaflag.composite
    
#    port = pd.read_csv('RESOLVE_full_snr5_dext_port.csv')#[portflag.sftoagn]
    port = pd.read_csv('ECO_raw_port.csv')#[portflag.sftoagn]
    port.index = port.name
    #jhu = pd.read_csv('RESOLVE_full_snr5_dext_jhu.csv')#[jhuflag.sftoagn]
    jhusel = pd.read_csv('ECO/SEL/ECO_full_snr5_dext_jhu.csv')#[jhuflag.sftoagn]
    jhusel.index = jhusel.name
    jhu = pd.read_csv('ECO_raw_jhu.csv')#[jhuflag.sftoagn]
#    jhu = pd.read_csv('RESOLVE_raw_jhu.csv')#[jhuflag.sftoagn]
    jhu.index = jhu.name
    nsasel = pd.read_csv('ECO/SEL/ECO_full_snr5_dext_nsa.csv')#.loc[nsaflag.index.values[nsaflag.sftoagn]]
    nsasel.index = nsasel.name
    nsa = pd.read_csv('ECO_raw_nsa.csv')
#    nsa = pd.read_csv('RESOLVE_raw_nsa.csv')
    nsa.index = nsa.name
    master = jhu#pd.read_csv('RESOLVE_snr5_master_bary.csv')
    master.index = master.name

if full: 
    eco = pd.read_csv('ECO+RESOLVE_inobssample.csv')
    eco.index = eco.name
    jhuflag = pd.read_csv('ECO/SEL/eco+resolve_emlineclass_dext_snr5_jhu.csv')
    jhuflag.index = jhuflag.galname
    jhuflagagn = jhuflag.sftoagn | jhuflag.agntosf | jhuflag.defagn | jhuflag.composite
    
    portflag = pd.read_csv('eco+resolve_emlineclass_dext_snr5_port.csv')
    portflag.index = portflag.galname
    portflagagn = portflag.sftoagn | portflag.agntosf | portflag.defagn | portflag.composite
    
    nsaflag = pd.read_csv('eco+resolve_emlineclass_dext_snr5_nsa.csv')
    nsaflag.index = nsaflag.galname
    nsaflagagn = nsaflag.sftoagn | nsaflag.agntosf | nsaflag.defagn | nsaflag.composite
    
    port = pd.read_csv('ECO+RESOLVE_snr5_dext_port.csv')#[portflag.sftoagn]
    port.index = port.name
    #jhu = pd.read_csv('RESOLVE_full_snr5_dext_jhu.csv')#[jhuflag.sftoagn]
    jhu = pd.read_csv('ECO+RESOLVE_snr5_dext_jhu.csv')#[jhuflag.sftoagn]
    jhu.index = jhu.name
    #nsa = pd.read_csv('RESOLVE_full_snr5_dext_nsa.csv')#.loc[nsaflag.index.values[nsaflag.sftoagn]]
    nsa = pd.read_csv('ECO+RESOLVE_snr5_dext_nsa.csv')
    nsa.index = nsa.name
    master = jhu#pd.read_csv('RESOLVE_snr5_master_bary.csv')
    master.index = master.name    


'''
plt.figure()
plt.plot(jhu.radeg, jhu.dedeg, 'bo', alpha = 0.3, mec = 'none')
plt.xlabel('RA')
plt.ylabel('Dec')
plt.figure()
plt.plot(jhu.radeg, jhu.dedeg, 'bo', alpha = 0.3, mec = 'none')
plt.plot(port.radeg, port.dedeg, 'go', alpha = 0.3, mec = 'none')
plt.xlabel('RA')
plt.ylabel('Dec')
plt.figure()
plt.plot(jhu.radeg, jhu.dedeg, 'bo', alpha = 0.3, mec = 'none')
plt.plot(port.radeg, port.dedeg, 'go', alpha = 0.3, mec = 'none')
plt.plot(nsa.radeg, nsa.dedeg, 'ro', alpha = 0.3, mec = 'none')
plt.xlabel('RA')
plt.ylabel('Dec')
plt.figure()
plt.plot(jhu.radeg, jhu.dedeg, 'bo', alpha = 0.3, mec = 'none')
plt.plot(nsa.radeg, nsa.dedeg, 'ro', alpha = 0.3, mec = 'none')
plt.xlabel('RA')
plt.ylabel('Dec')
'''
ecoinjhu = np.array([x in jhu.name for x in eco.name])
ecoinport = np.array([x in port.name for x in eco.name])
ecoinnsa = np.array([x in nsa.name for x in eco.name])

port_jhu = (ecoinjhu & ecoinport)
jhu_nsa = (ecoinjhu & ecoinnsa)
port_nsa = (ecoinnsa & ecoinport)

jhu_nsanames = jhu.name.loc[eco.name[jhu_nsa]]
jhuagn = jhuflag.galname[~jhuflag.defstarform]
nsaagn = nsaflag.galname[~nsaflag.defstarform]
jhunsasel = np.intersect1d(nsasel.name, jhusel.name)
jhunsasel = np.intersect1d(jhunsasel, jhu_nsanames)
jhunsaagn = np.intersect1d(nsaagn, jhuagn)
jhunsaagn = np.intersect1d(jhunsaagn, jhu_nsanames)

print('Total Number of SEL ECO galaxies in SDSS:'+str(len(master)))
print('Number of SEL ECO galaxies in each catalog:')
print('\nJHU:'+str(len(jhu)))
print('\nPort:'+str(len(port)))
print('\nNSA:'+str(len(nsa)))

print('\n\nOverlap between catalogs: \n')
print('\nJHU and Port:'+str(np.sum(jhu_port)))
print('\nJHU and NSA:'+str(np.sum(jhu_nsa)))
print('\nNSA and Port:'+str(np.sum(nsa_port)))

lines = ['h_alpha_flux','h_beta_flux', 'nii_6584_flux', 'oi_6300_flux', 
         'oiii_5007_flux', 'sii_6731_flux']
#lines = ['nii_6584_flux']
#lines = ['h_alpha_flux_err','h_beta_flux_err', 'nii_6584_flux_err', 'oi_6300_flux_err', 
#         'oiii_5007_flux_err', 'sii_6731_flux_err']
errlines = ['h_alpha_flux_err','h_beta_flux_err', 'nii_6584_flux_err', 'oi_6300_flux_err', 
         'oiii_5007_flux_err', 'sii_6731_flux_err']
numlines = ['nii_6584_flux', 'oi_6300_flux', 
         'oiii_5007_flux', 'sii_6731_flux']
denlines = ['h_alpha_flux','h_alpha_flux','h_beta_flux','h_alpha_flux'] 

portlines = {'h_alpha_flux':'Flux_Ha_6562',
         'h_beta_flux': 'Flux_Hb_4861', 
         'nii_6584_flux': 'Flux_NII_6583', 
         'oi_6300_flux': 'Flux_OI_6300', 
         'oiii_5007_flux': 'Flux_OIII_5006', 
         'sii_6731_flux': 'Flux_SII_6730'}

blue = (jhu.logmstar < 9.5) & (jhu.modelu_rcorr < 1.5)
fig, axes = plt.subplots(2,3, figsize = (13.5,6.5))
fig.suptitle("NSA vs. JHU")
for line in lines:
    ax = axes.ravel()[lines.index(line)]
    if line == 'sii_6731_flux':
        jhuline = jhu[line].loc[eco.name[jhu_nsa]] + jhu['sii_6717_flux'].loc[eco.name[jhu_nsa]]
#        jhuline = np.sqrt(jhu[line].loc[eco.name[jhu_nsa]]**2 + jhu['sii_6717_flux_err'].loc[eco.name[jhu_nsa]]**2)
    else:
        jhuline = jhu[line].loc[eco.name[jhu_nsa]] 
    nsaline = nsa[line].loc[eco.name[jhu_nsa]]
    good = (np.isfinite(jhuline) & (jhuline> 0) & (jhuline< 1e5) & \
           np.isfinite(nsaline) & (nsaline> 0) & (nsaline< 1e5))
    xline = jhuline[good]
    yline = nsaline[good]
    xvalue = 10**np.log10(xline)
    yvalue = 10**np.log10(yline)
    residual = (yvalue - xvalue)/xvalue
    ax.plot(xvalue, residual, 'k.')
    ax.plot(xvalue[good & blue], residual[good & blue], 'b.')
    ax.set_xlabel('JHU '+line)
    ax.set_ylabel('(NSA - JHU)/JHU')
    xaxis = np.arange(np.min(xvalue),np.max(xvalue),0.1)
    ax.plot(xaxis, 0*xaxis,'r')
    medoffset = np.nanmedian(residual[np.isfinite(residual)])
    ax.plot(xaxis, 0*xaxis+medoffset,'r-.')
    ax.set_ylim(-0.5,5)
    ax.set_xscale('log')
#manager = plt.get_current_fig_manager()
#manager.window.showMaximized()
fig.tight_layout()#w_pad = 0.1)

fig, axes = plt.subplots(2,3, figsize = (13.5,6.5))
for line in lines:
    ax = axes.ravel()[lines.index(line)]
    if line == 'sii_6731_flux':
        jhuline = jhu[line+'_err'].loc[eco.name[jhu_nsa]] + jhu['sii_6717_flux'].loc[eco.name[jhu_nsa]]
        jhulinerr = np.sqrt(jhu[line+'_err'].loc[eco.name[jhu_nsa]]**2 + jhu['sii_6717_flux_err'].loc[eco.name[jhu_nsa]]**2)
    else:
        jhuline = jhu[line].loc[eco.name[jhu_nsa]] 
        jhulineerr = jhu[line+'_err'].loc[eco.name[jhu_nsa]] 
    nsaline = nsa[line].loc[eco.name[jhu_nsa]]
    nsalineerr = nsa[line+'_err'].loc[eco.name[jhu_nsa]]
    good = (np.isfinite(jhuline) & (jhuline> 0) & (jhuline< 1e2) & \
           np.isfinite(nsaline) & (nsaline> 0) & (nsaline< 1e2))

    xvalue = jhulineerr[good]
    yvalue = nsalineerr[good]

    residual = (yvalue/xvalue) #- xvalue)/xvalue

#    good = residual>
    mass = jhu.logmstar.loc[good.index[good]]
#    mass = mass[residual<10]
#    residual = residual[residual<10]
    ax.plot(mass, residual, 'k.')
    ax.plot(mass[good & blue], residual[good & blue], 'b.')
    goodsel = np.intersect1d(jhunsasel, eco.name[jhu_nsa][good])
    ax.plot(mass.loc[goodsel], residual.loc[goodsel], 'r.')
    goodagn = np.intersect1d(jhunsaagn, eco.name[jhu_nsa][good])
    ax.plot(mass.loc[goodagn], residual.loc[goodagn], 'c.')
    ax.set_title('JHU '+line)
    ax.set_ylabel('(NSA/JHU)')
    ax.set_xlabel('stellar mass')
    xaxis = np.arange(np.min(mass),np.max(mass),0.1)
    ax.plot(xaxis, 0*xaxis,'r')
    medoffset = np.nanmedian(residual[np.isfinite(residual)])
    ax.plot(xaxis, 0*xaxis+medoffset,'r-.')
    ax.set_ylim(-0.5,5)
    print(line+'_err', stats.spearmanr(10**mass, residual))
    cc, pnull = stats.spearmanr(10**mass, residual)
#manager = plt.get_current_fig_manager()
#manager.window.showMaximized()
fig.tight_layout()#w_pad = 0.1)

fig, axes = plt.subplots(2,3, figsize = (13.5,6.5))
for line in lines:
    ax = axes.ravel()[lines.index(line)]
    if line == 'sii_6731_flux':
        jhuline = jhu[line+'_err'].loc[eco.name[jhu_nsa]] + jhu['sii_6717_flux'].loc[eco.name[jhu_nsa]]
        jhulinerr = np.sqrt(jhu[line+'_err'].loc[eco.name[jhu_nsa]]**2 + jhu['sii_6717_flux_err'].loc[eco.name[jhu_nsa]]**2)
    else:
        jhuline = jhu[line].loc[eco.name[jhu_nsa]] 
        jhulineerr = jhu[line+'_err'].loc[eco.name[jhu_nsa]] 
    nsaline = nsa[line].loc[eco.name[jhu_nsa]]
    nsalineerr = nsa[line+'_err'].loc[eco.name[jhu_nsa]]
    good = (np.isfinite(jhuline) & (jhuline> 0) & (jhuline< 1e5) & \
           np.isfinite(nsaline) & (nsaline> 0) & (nsaline< 1e5))
    xline = jhuline[good]
    yline = nsaline[good]
    xvalue = xline/jhulineerr[good]
    yvalue = yline/nsalineerr[good]

    residual = (yvalue - xvalue)/xvalue
    mass = jhu.logmstar.loc[good.index[good]]
    ax.plot(mass, residual, 'k.')
    ax.plot(mass[good & blue], residual[good & blue], 'b.')
    goodsel = np.intersect1d(jhunsasel, eco.name[jhu_nsa][good])
    ax.plot(mass.loc[goodsel], residual.loc[goodsel], 'r.')
    goodagn = np.intersect1d(jhunsaagn, eco.name[jhu_nsa][good])
    ax.plot(mass.loc[goodagn], residual.loc[goodagn], 'c.')
    ax.set_title('JHU '+line)
    ax.set_ylabel('(NSA - JHU)/JHU')
    ax.set_xlabel('stellar mass')
    xaxis = np.arange(np.min(mass),np.max(mass),0.1)
    ax.plot(xaxis, 0*xaxis,'r')
    medoffset = np.nanmedian(residual[np.isfinite(residual)])
    ax.plot(xaxis, 0*xaxis+medoffset,'r-.')
    ax.set_ylim(-0.5,5)
#    print(line, stats.spearmanr(mass, residual))

#manager = plt.get_current_fig_manager()
#manager.window.showMaximized()
fig.tight_layout()#w_pad = 0.1)

fig, axes = plt.subplots(2,3, figsize = (13.5,6.5))
for line in lines:
    ax = axes.ravel()[lines.index(line)]
    portline = port[line].loc[eco.name[port_nsa]] 
    portlineerr = port[line+'_err'].loc[eco.name[port_nsa]] 
    nsaline = nsa[line].loc[eco.name[port_nsa]]
    nsalineerr = nsa[line+'_err'].loc[eco.name[port_nsa]]
    good = (np.isfinite(portline) & (portline> 0) & (portline< 1e5) & \
           np.isfinite(nsaline) & (nsaline> 0) & (nsaline< 1e5))
    xline = portline[good]
    yline = nsaline[good]
#    xvalue = xline/portlineerr[good]
#    yvalue = yline/nsalineerr[good]

    xvalue = portlineerr[good]
    yvalue = nsalineerr[good]
    residual = (yvalue - xvalue)/xvalue
    mass = port.logmstar.loc[good.index[good]]
    ax.plot(mass, residual, 'k.')
    ax.plot(mass[good & blue], residual[good & blue], 'b.')
    ax.set_title('port '+line)
    ax.set_ylabel('(NSA - port)/port')
    ax.set_xlabel('stellar mass')
    xaxis = np.arange(np.min(mass),np.max(mass),0.1)
    ax.plot(xaxis, 0*xaxis,'r')
    medoffset = np.nanmedian(residual[np.isfinite(residual)])
    ax.plot(xaxis, 0*xaxis+medoffset,'r-.')
    ax.set_ylim(-0.5,5)
fig.tight_layout()

fig, axes = plt.subplots(2,3, figsize = (13.5,6.5))
for line in lines:
    ax = axes.ravel()[lines.index(line)]
    portline = port[portlines[line]].loc[eco.name[port_nsa]] 
    portlineerr = port[portlines[line]+'_Err'].loc[eco.name[port_nsa]] 
    jhuline = jhu[line].loc[eco.name[port_jhu]]
    jhulineerr = jhu[line+'_err'].loc[eco.name[port_jhu]]
    good = (np.isfinite(portline) & (portline> 0) & (portline< 1e5) & \
           np.isfinite(jhuline) & (jhuline> 0) & (jhuline< 1e5))
    xline = portline[good]
    yline = jhuline[good]
#    xvalue = xline/portlineerr[good]
#    yvalue = yline/jhulineerr[good]

    yvalue = portlineerr[good]
    xvalue = jhulineerr[good]
    residual = yvalue/xvalue
    mass = port.logmstar.loc[good.index[good]]
    ax.plot(mass, residual, 'k.')
    ax.plot(mass[good & blue], residual[good & blue], 'b.')
    ax.set_title('jhu '+line)
    ax.set_ylabel('(port - jhu)/jhu')
    ax.set_xlabel('stellar mass')
    xaxis = np.arange(np.min(mass),np.max(mass),0.1)
    ax.plot(xaxis, 0*xaxis,'r')
    medoffset = np.nanmedian(residual[np.isfinite(residual)])
    ax.plot(xaxis, 0*xaxis+medoffset,'r-.')
    ax.set_ylim(-0.5,5)
fig.tight_layout()

#fig, axes = plt.subplots(2,3, figsize = (13.5,6.5))
#fig.suptitle("NSA vs. port")
#for line in lines:
#    ax = axes.ravel()[lines.index(line)]
#    if line == 'sii_6731_flux':
#        xvalue = jhu[line].loc[eco.name[jhu_nsa]] + jhu['sii_6717_flux'].loc[eco.name[jhu_nsa]]
#    else:
#        xvalue = jhu[line].loc[eco.name[jhu_nsa]] 
#    yvalue = 10**np.log10(nsa[line].loc[eco.name[jhu_nsa]])
#    residual = (yvalue - xvalue)/xvalue
#    mass = jhu.logmstar.loc[eco.name[jhu_nsa]]
#    ax.plot(mass, residual, 'b.')
#    ax.set_xlabel('JHU '+line)
#    ax.set_ylabel('(NSA - JHU)/JHU')
#    xaxis = np.arange(np.min(mass),np.max(mass),0.1)
#    ax.plot(xaxis, 0*xaxis,'k')
#    medoffset = np.nanmedian(residual[np.isfinite(residual)])
#    ax.plot(xaxis, 0*xaxis+medoffset,'k-.')
#    ax.set_ylim(-0.5,5)
##manager = plt.get_current_fig_manager()
##manager.window.showMaximized()
#fig.tight_layout()#w_pad = 0.1)

fig, axes = plt.subplots(2,2, figsize = (13.5,6.5))
fig.suptitle("NSA vs. JHU")
for i in range(len(numlines)):
    
    numline = numlines[i]
    denline = denlines[i]
    ax = axes.ravel()[numlines.index(numline)]
    
    if i == 3: #numline == 'sii_6731_flux':
        print( 'added')
        jhunum = jhu[numline].loc[eco.name[jhu_nsa]]+jhu['sii_6717_flux'].loc[eco.name[jhu_nsa]]
    else:
        jhunum = jhu[numline].loc[eco.name[jhu_nsa]]
 
    jhuden = jhu[denline].loc[eco.name[jhu_nsa]]
    xvalue = 10**np.log10(jhunum/jhuden)
    print(numline)
    nsanum = nsa[numline].loc[eco.name[jhu_nsa]]
    nsaden = nsa[denline].loc[eco.name[jhu_nsa]]
    yvalue = 10**np.log10(nsanum/nsaden)
    residual = (yvalue - xvalue)/xvalue
    ax.plot(xvalue, residual, 'r.')
    ax.set_xlabel('JHU '+numline+'/'+denline)
    ax.set_ylabel('(NSA - JHU)/JHU')#+numline+'/'+denline)
    xaxis = np.arange(np.min(xvalue[np.isfinite(xvalue)]),np.max(xvalue[np.isfinite(xvalue)]),0.1)
    ax.plot(xaxis, 0*xaxis,'k')
    medoffset = np.nanmedian(residual[np.isfinite(residual)])
    ax.plot(xaxis, 0*xaxis+medoffset,'k-.')
    print(np.nanmedian(residual[np.isfinite(residual)]))
    print(min(residual), max(residual))
    ax.set_ylim(-0.5,+2)
#manager = plt.get_current_fig_manager()
#manager.window.showMaximized()
fig.tight_layout()

fig, axes = plt.subplots(2,2, figsize = (13.5,6.5))
fig.suptitle("NSA vs. JHU")
for i in range(len(numlines)):
    
    numline = numlines[i]
    denline = denlines[i]
    ax = axes.ravel()[numlines.index(numline)]
    jhunum = jhu[numline].loc[eco.name[jhu_nsa]]
    jhuden = jhu[denline].loc[eco.name[jhu_nsa]]
    xvalue = 10**np.log10(jhunum/jhuden)
    print(numline)
    nsanum = nsa[numline].loc[eco.name[jhu_nsa]]
    nsaden = nsa[numline].loc[eco.name[jhu_nsa]]
    yvalue = 10**np.log10(nsanum/nsaden)
    residual = (yvalue - xvalue)/xvalue
    ax.plot(xvalue, residual, 'r.')
    ax.set_xlabel('JHU '+numline+'/'+denline)
    ax.set_ylabel('(NSA - JHU)/JHU')#+numline+'/'+denline)
    xaxis = np.arange(np.min(xvalue[np.isfinite(xvalue)]),np.max(xvalue[np.isfinite(xvalue)]),0.1)
    ax.plot(xaxis, 0*xaxis,'k')
    print(np.nanmedian(residual[np.isfinite(residual)]))
    print(min(residual), max(residual))
    #ax.set_ylim(-1,+2)
#manager = plt.get_current_fig_manager()
#manager.window.showMaximized()
fig.tight_layout()

#fig, axes = plt.subplots(2,3, figsize = (13.5,6.5))
#fig.suptitle("Portsmouth vs. JHU")
#for line in lines:
#    ax = axes.ravel()[lines.index(line)]
#    xvalue = jhu[line].loc[eco.name[jhu_port]]
#    yvalue = port[portlines[line]].loc[eco.name[jhu_port]]
#    residual = (yvalue - xvalue)/xvalue
#    ax.plot(xvalue, residual, 'r.')
#    ax.set_xlabel('JHU '+line)
#    ax.set_ylabel('Portsmouth '+line)
#    xaxis = np.arange(np.max(jhu[line].loc[eco.name[jhu_nsa]]))
#    ax.plot(xaxis, 0*xaxis,'k')
##manager = plt.get_current_fig_manager()
##manager.window.showMaximized()
#fig.tight_layout()

#fig, axes = plt.subplots(2,2, figsize = (13.5,6.5))
#fig.suptitle("Portsmouth vs. JHU")
#for i in range(len(numlines)):
#    numline = numlines[i]
#    denline = denlines[i]
#    ax = axes.ravel()[numlines.index(numline)]
#    xvalue = np.log10(jhu[numline].loc[eco.name[jhu_port]]/jhu[denline].loc[eco.name[jhu_port]])
#    yvalue = np.log10(port[portlines[numline]].loc[eco.name[jhu_port]]/port[portlines[denline]].loc[eco.name[jhu_port]])
#    residual = (yvalue - xvalue)/xvalue
#    ax.plot(xvalue, yvalue, 'r.')
#    ax.set_xlabel('JHU '+line)
#    ax.set_ylabel('Portsmouth '+line)
#    xaxis = np.arange(np.min(xvalue),np.max(xvalue),0.1)
#    ax.plot(xaxis, xaxis,'k')
##manager = plt.get_current_fig_manager()
##manager.window.showMaximized()
#fig.tight_layout()
#
#
#fig, axes = plt.subplots(2,3, figsize = (13.5,6.5))
#fig.suptitle("NSA vs. Portsmouth")
#for line in lines:
#    ax = axes.ravel()[lines.index(line)]
#    ax.plot(port[portlines[line]].loc[eco.name[nsa_port]], 
#            nsa[line].loc[eco.name[nsa_port]], 'g.')
#    ax.set_xlabel('Portsmouth '+line)
#    ax.set_ylabel('NSA '+line)
#    xaxis = np.arange(np.max(nsa[line].loc[eco.name[jhu_nsa]]))
#    ax.plot(xaxis, xaxis,'k')
##manager = plt.get_current_fig_manager()
##manager.window.showMaximized()
#fig.tight_layout()
#
#conf = pd.read_csv('ECO_snr5_master_conf.csv')
#conf.index = conf.name
#conf1 = conf.name[conf.confidence_level == 1.0]
#print("ECO galaxies with AGN confidence level = 1.0 :"+ str(len(conf1)))
#conf1_nsa = master.loc[conf1][master.loc[conf1].source == 'nsa']
#conf1_jhu = master.loc[conf1][master.loc[conf1].source == 'jhu']
#conf1_port = master.loc[conf1][master.loc[conf1].source == 'port']
#
#print("\n\n Confidence level 1 AGN from JHU :"+ str(len(conf1_jhu)))
#print("\n\n Confidence level 1 AGN from Portsmouth :"+ str(len(conf1_port)))
#print("\n\n Confidence level 1 AGN from NSA :"+ str(len(conf1_nsa)))


ratiolines = ['nii_6584_flux', 'oi_6300_flux', 
         'oiii_5007_flux', 'sii_6731_flux']

fig, axes = plt.subplots(2,2, figsize = (13.5,6.5))
fig.suptitle("NSA vs. JHU Residuals")
for line in ratiolines:
    if line != ratiolines[2]:
        x_ratio = jhu[line].loc[eco.name[jhu_nsa]]/jhu[lines[0]].loc[eco.name[jhu_nsa]]
        y_ratio = nsa[line].loc[eco.name[jhu_nsa]]/nsa[lines[0]].loc[eco.name[jhu_nsa]]
        residual = (x_ratio - y_ratio)/x_ratio
        ax = axes.ravel()[ratiolines.index(line)]
        ax.plot(x_ratio,residual, 'b.')
        ax.set_xlabel('JHU '+line)
        ax.set_ylabel('NSA '+line)
        xaxis = np.arange(0, np.max(jhu[line].loc[eco.name[jhu_nsa]]/jhu[lines[0]].loc[eco.name[jhu_nsa]]),0.01)
        ax.plot(xaxis, np.zeros(len(xaxis)),'k')
#manager = plt.get_current_fig_manager()
#manager.window.showMaximized()
fig.tight_layout()#w_pad = 0.1)


#reineslines = {'h_alpha_flux':'Ha-n',
#         'h_beta_flux': 'Hb-n', 
#         'nii_6584_flux': 'NII-6583', 
#         'oi_6300_flux': 'OI-6300', 
#         'oiii_5007_flux': 'OIII-5007', 
#         'sii_6731_flux': 'SII-6731',
#         'sii_6716_flux': 'SII-6716'}
#
#nsalines = {'h_alpha_flux':'HAFLUX',
#         'h_beta_flux': 'HBFLUX', 
#         'nii_6584_flux': 'N2FLUX', 
#         'oi_6300_flux': 'O1FLUX', 
#         'oiii_5007_flux': 'O3FLUX', 
#         'sii_6731_flux': 'S2FLUX'}
#
#reines = pd.read_csv("Reines13flux.csv")
#reines.index = reines.NSAID
#from astropy.table import Table
#nsa = Table.read('nsa_v0_1_2.fits')
#
#ndx = [x for x in range(len(nsa)) if nsa['NSAID'][x] in reines.NSAID]
#
#fig, axes = plt.subplots(2,3, figsize = (13.5,6.5))
#fig.suptitle("NSA vs. Reines")
#for line in lines:
#    ax = axes.ravel()[lines.index(line)]
#    if line == 'sii_6731_flux':        
#        xvalue = 10**np.log10(reines[reineslines[line]].loc[nsa['NSAID'][ndx]]+\
#                          reines['SII-6716'].loc[nsa['NSAID'][ndx]]) 
#    else:
#        xvalue = 10**np.log10(reines[reineslines[line]].loc[nsa['NSAID'][ndx]])
#        
#    yvalue = 10**np.log10(nsa[nsalines[line]][ndx])
#    residual = (yvalue - xvalue)/xvalue
#    #xvalue = np.log10(xvalue)
#    ax.plot(xvalue, residual,'g.')
#    ax.set_xlabel('Reines '+line)
#    ax.set_ylabel('(NSA - Reines)/Reines')
#    ax.set_ylim(-0.2,1)
#    xaxis = np.arange(np.min(xvalue[np.isfinite(xvalue)]),np.max(xvalue[np.isfinite(xvalue)]),0.1)
#    ax.set_xscale('log')
#    ax.plot(xaxis, 0*xaxis,'k')
#    medoffset = np.nanmedian(residual[np.isfinite(residual)])
#    ax.plot(xaxis, 0*xaxis+medoffset,'k-.')
#    
##manager = plt.get_current_fig_manager()
##manager.window.showMaximized()
#fig.tight_layout()
#
#fig, axes = plt.subplots(2,2, figsize = (13.5,6.5))
#fig.suptitle("NSA vs. Reines")
#for i in range(len(numlines)):
#    
#    numline = numlines[i]
#    denline = denlines[i]
#    ax = axes.ravel()[numlines.index(numline)]
#    if numline == 'sii_6731_flux':        
#        reinesnum = (reines[reineslines[line]].loc[nsa['NSAID'][ndx]]+\
#                          reines['SII-6716'].loc[nsa['NSAID'][ndx]]) 
#    else:
#        reinesnum = reines[reineslines[numline]].loc[nsa['NSAID'][ndx]]
#    reinesden = reines[reineslines[denline]].loc[nsa['NSAID'][ndx]]
#    xvalue = 10**np.log10(reinesnum/reinesden)
#    print(numline)
#    nsanum = nsa[nsalines[numline]][ndx]
#    nsaden = nsa[nsalines[denline]][ndx]
#    yvalue = 10**np.log10(nsanum/nsaden)
#    residual = (yvalue - xvalue)/xvalue
#    ax.plot(xvalue, residual, 'r.')
#    ax.set_xlabel('Reines '+numline+'/'+denline)
#    ax.set_ylabel('(NSA - Reines)/Reines')#+numline+'/'+denline)
#    ax.set_ylim(-0.2,1)
#    ax.set_xscale('log')
#    xaxis = np.arange(np.min(xvalue[np.isfinite(xvalue)]),np.max(xvalue[np.isfinite(xvalue)]),0.1)
#    ax.plot(xaxis, 0*xaxis,'k')
#    medoffset = np.nanmedian(residual[np.isfinite(residual)])
#    ax.plot(xaxis, 0*xaxis+medoffset,'k-.')
#    print(np.nanmedian(residual[np.isfinite(residual)]))
#    print(min(residual), max(residual))
##manager = plt.get_current_fig_manager()
##manager.window.showMaximized()
#fig.tight_layout()
