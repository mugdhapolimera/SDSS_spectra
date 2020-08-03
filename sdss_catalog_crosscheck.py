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

#eco = pd.read_csv('ECO_live22Oct2018.csv')
#eco.index = eco.name
#jhuflag = pd.read_csv('eco_emlineclass_full_bary_jhu.csv')
#jhuflag.index = jhuflag.galname
#portflag = pd.read_csv('eco_emlineclass_full_snr5_port.csv')
#portflag.index = portflag.galname
#nsaflag = pd.read_csv('eco_emlineclass_full_snr5_nsa.csv')
#nsaflag.index = nsaflag.galname
#
#port = pd.read_csv('ECO_full_bary_port.csv')#[portflag.sftoagn]
#port.index = port.name
#jhu = pd.read_csv('ECO_full_bary_jhu.csv')#[jhuflag.sftoagn]
#jhu.index = jhu.name
#nsa = pd.read_csv('NSA_ECO_snr5.csv')#.loc[nsaflag.index.values[nsaflag.sftoagn]]
#nsa.index = nsa.name
#master = pd.read_csv('ECO_snr5_master_bary.csv')
#master.index = master.name

eco = pd.read_csv('RESOLVE_liveMay2020.csv')
eco.index = eco.name
jhuflag = pd.read_csv('resolve_emlineclass_dext_snr5_jhu.csv')
jhuflag.index = jhuflag.galname
portflag = pd.read_csv('resolve_emlineclass_dext_snr5_port.csv')
portflag.index = portflag.galname
nsaflag = pd.read_csv('resolve_emlineclass_dext_snr5_nsa.csv')
nsaflag.index = nsaflag.galname

port = pd.read_csv('RESOLVE_full_snr5_dext_port.csv')#[portflag.sftoagn]
port.index = port.name
#jhu = pd.read_csv('RESOLVE_full_snr5_dext_jhu.csv')#[jhuflag.sftoagn]
jhu = pd.read_csv('RESOLVE_full_snr5_dext_jhu.csv')#[jhuflag.sftoagn]
jhu.index = jhu.name
#nsa = pd.read_csv('RESOLVE_full_snr5_dext_nsa.csv')#.loc[nsaflag.index.values[nsaflag.sftoagn]]
nsa = pd.read_csv('RESOLVE_full_snr5_dext_nsa.csv')
nsa.index = nsa.name
master = jhu#pd.read_csv('RESOLVE_snr5_master_bary.csv')
#master.index = master.name
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

jhu_port = (ecoinjhu & ecoinport)
jhu_nsa = (ecoinjhu & ecoinnsa)
nsa_port = (ecoinnsa & ecoinport)

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
numlines = ['nii_6584_flux', 'oi_6300_flux', 
         'oiii_5007_flux', 'sii_6731_flux']
denlines = ['h_alpha_flux','h_alpha_flux','h_beta_flux','h_alpha_flux'] 

portlines = {'h_alpha_flux':'Flux_Ha_6562',
         'h_beta_flux': 'Flux_Hb_4861', 
         'nii_6584_flux': 'Flux_NII_6583', 
         'oi_6300_flux': 'Flux_OI_6300', 
         'oiii_5007_flux': 'Flux_OIII_5006', 
         'sii_6731_flux': 'Flux_SII_6730'}

fig, axes = plt.subplots(2,3, figsize = (13.5,6.5))
fig.suptitle("NSA vs. JHU")
for line in lines:
    ax = axes.ravel()[lines.index(line)]
    xvalue = 10**np.log10(jhu[line].loc[eco.name[jhu_nsa]])
    yvalue = 10**np.log10(nsa[line].loc[eco.name[jhu_nsa]])
    residual = (yvalue - xvalue)/xvalue
    ax.plot(xvalue, residual, 'b.')
    ax.set_xlabel('JHU '+line)
    ax.set_ylabel('(NSA - JHU)/JHU')
    xaxis = np.arange(np.min(xvalue),np.max(xvalue),0.1)
    ax.plot(xaxis, 0*xaxis,'k')
    medoffset = np.nanmedian(residual[np.isfinite(residual)])
    ax.plot(xaxis, 0*xaxis+medoffset,'k-.')
    ax.set_ylim(-0.5,5)
#manager = plt.get_current_fig_manager()
#manager.window.showMaximized()
fig.tight_layout()#w_pad = 0.1)

fig, axes = plt.subplots(2,3, figsize = (13.5,6.5))
fig.suptitle("NSA vs. JHU")
for line in lines:
    ax = axes.ravel()[lines.index(line)]
    if line == 'sii_6731_flux':
        xvalue = jhu[line].loc[eco.name[jhu_nsa]] + jhu['sii_6717_flux'].loc[eco.name[jhu_nsa]]
    else:
        xvalue = jhu[line].loc[eco.name[jhu_nsa]] 
    yvalue = 10**np.log10(nsa[line].loc[eco.name[jhu_nsa]])
    residual = (yvalue - xvalue)/xvalue
    ax.plot(xvalue, residual*100.0,'b.')
    ax.set_xlabel('JHU '+line)
    ax.set_ylabel('NSA '+line)
    xaxis = np.arange(np.min(xvalue),np.max(xvalue),0.1)
    ax.plot(xaxis, 0*xaxis,'k')
    medoffset = np.nanmedian(residual[np.isfinite(residual)])
    #ax.plot(xaxis, 0*xaxis+medoffset,'k-.')
    ax.set_ylim(-0.5,200)
#manager = plt.get_current_fig_manager()
#manager.window.showMaximized()
fig.tight_layout()#w_pad = 0.1)

fig, axes = plt.subplots(2,3, figsize = (13.5,6.5))
fig.suptitle("NSA vs. JHU")
for line in lines:
    ax = axes.ravel()[lines.index(line)]
    if line == 'sii_6731_flux':
        xvalue = jhu[line].loc[eco.name[jhu_nsa]] + jhu['sii_6717_flux'].loc[eco.name[jhu_nsa]]
    else:
        xvalue = jhu[line].loc[eco.name[jhu_nsa]] 
    yvalue = 10**np.log10(nsa[line].loc[eco.name[jhu_nsa]])
    residual = (yvalue - xvalue)/xvalue
    mass = jhu.logmstar.loc[eco.name[jhu_nsa]]
    ax.plot(mass, residual, 'b.')
    ax.set_xlabel('JHU '+line)
    ax.set_ylabel('(NSA - JHU)/JHU')
    xaxis = np.arange(np.min(mass),np.max(mass),0.1)
    ax.plot(xaxis, 0*xaxis,'k')
    medoffset = np.nanmedian(residual[np.isfinite(residual)])
    ax.plot(xaxis, 0*xaxis+medoffset,'k-.')
    ax.set_ylim(-0.5,5)
#manager = plt.get_current_fig_manager()
#manager.window.showMaximized()
fig.tight_layout()#w_pad = 0.1)

fig, axes = plt.subplots(2,2, figsize = (13.5,6.5))
fig.suptitle("NSA vs. JHU")
for i in range(len(numlines)):
    
    numline = numlines[i]
    denline = denlines[i]
    ax = axes.ravel()[numlines.index(numline)]
    
    if i == 3: #numline == 'sii_6731_flux':
        print 'added'
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

fig, axes = plt.subplots(2,3, figsize = (13.5,6.5))
fig.suptitle("Portsmouth vs. JHU")
for line in lines:
    ax = axes.ravel()[lines.index(line)]
    xvalue = jhu[line].loc[eco.name[jhu_port]]
    yvalue = port[portlines[line]].loc[eco.name[jhu_port]]
    residual = (yvalue - xvalue)/xvalue
    ax.plot(xvalue, residual, 'r.')
    ax.set_xlabel('JHU '+line)
    ax.set_ylabel('Portsmouth '+line)
    xaxis = np.arange(np.max(jhu[line].loc[eco.name[jhu_nsa]]))
    ax.plot(xaxis, 0*xaxis,'k')
#manager = plt.get_current_fig_manager()
#manager.window.showMaximized()
fig.tight_layout()

fig, axes = plt.subplots(2,2, figsize = (13.5,6.5))
fig.suptitle("Portsmouth vs. JHU")
for i in range(len(numlines)):
    numline = numlines[i]
    denline = denlines[i]
    ax = axes.ravel()[numlines.index(numline)]
    xvalue = np.log10(jhu[numline].loc[eco.name[jhu_port]]/jhu[denline].loc[eco.name[jhu_port]])
    yvalue = np.log10(port[portlines[numline]].loc[eco.name[jhu_port]]/port[portlines[denline]].loc[eco.name[jhu_port]])
    residual = (yvalue - xvalue)/xvalue
    ax.plot(xvalue, yvalue, 'r.')
    ax.set_xlabel('JHU '+line)
    ax.set_ylabel('Portsmouth '+line)
    xaxis = np.arange(np.min(xvalue),np.max(xvalue),0.1)
    ax.plot(xaxis, xaxis,'k')
#manager = plt.get_current_fig_manager()
#manager.window.showMaximized()
fig.tight_layout()


fig, axes = plt.subplots(2,3, figsize = (13.5,6.5))
fig.suptitle("NSA vs. Portsmouth")
for line in lines:
    ax = axes.ravel()[lines.index(line)]
    ax.plot(port[portlines[line]].loc[eco.name[nsa_port]], 
            nsa[line].loc[eco.name[nsa_port]], 'g.')
    ax.set_xlabel('Portsmouth '+line)
    ax.set_ylabel('NSA '+line)
    xaxis = np.arange(np.max(nsa[line].loc[eco.name[jhu_nsa]]))
    ax.plot(xaxis, xaxis,'k')
#manager = plt.get_current_fig_manager()
#manager.window.showMaximized()
fig.tight_layout()

conf = pd.read_csv('ECO_snr5_master_conf.csv')
conf.index = conf.name
conf1 = conf.name[conf.confidence_level == 1.0]
print("ECO galaxies with AGN confidence level = 1.0 :"+ str(len(conf1)))
conf1_nsa = master.loc[conf1][master.loc[conf1].source == 'nsa']
conf1_jhu = master.loc[conf1][master.loc[conf1].source == 'jhu']
conf1_port = master.loc[conf1][master.loc[conf1].source == 'port']

print("\n\n Confidence level 1 AGN from JHU :"+ str(len(conf1_jhu)))
print("\n\n Confidence level 1 AGN from Portsmouth :"+ str(len(conf1_port)))
print("\n\n Confidence level 1 AGN from NSA :"+ str(len(conf1_nsa)))

from matplotlib_venn import venn3
set1 = set(jhu.name)
set2 = set(nsa.name)
set3 = set(port.name)
plt.figure()
venn3([set1, set2, set3], ('JHU', 'NSA', 'Port'))

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


reineslines = {'h_alpha_flux':'Ha-n',
         'h_beta_flux': 'Hb-n', 
         'nii_6584_flux': 'NII-6583', 
         'oi_6300_flux': 'OI-6300', 
         'oiii_5007_flux': 'OIII-5007', 
         'sii_6731_flux': 'SII-6731',
         'sii_6716_flux': 'SII-6716'}

nsalines = {'h_alpha_flux':'HAFLUX',
         'h_beta_flux': 'HBFLUX', 
         'nii_6584_flux': 'N2FLUX', 
         'oi_6300_flux': 'O1FLUX', 
         'oiii_5007_flux': 'O3FLUX', 
         'sii_6731_flux': 'S2FLUX'}

reines = pd.read_csv("Reines13flux.csv")
reines.index = reines.NSAID
from astropy.table import Table
nsa = Table.read('nsa_v0_1_2.fits')

ndx = [x for x in range(len(nsa)) if nsa['NSAID'][x] in reines.NSAID]

fig, axes = plt.subplots(2,3, figsize = (13.5,6.5))
fig.suptitle("NSA vs. Reines")
for line in lines:
    ax = axes.ravel()[lines.index(line)]
    if line == 'sii_6731_flux':        
        xvalue = 10**np.log10(reines[reineslines[line]].loc[nsa['NSAID'][ndx]]+\
                          reines['SII-6716'].loc[nsa['NSAID'][ndx]]) 
    else:
        xvalue = 10**np.log10(reines[reineslines[line]].loc[nsa['NSAID'][ndx]])
        
    yvalue = 10**np.log10(nsa[nsalines[line]][ndx])
    residual = (yvalue - xvalue)/xvalue
    #xvalue = np.log10(xvalue)
    ax.plot(xvalue, residual,'g.')
    ax.set_xlabel('Reines '+line)
    ax.set_ylabel('(NSA - Reines)/Reines')
    ax.set_ylim(-0.2,1)
    xaxis = np.arange(np.min(xvalue[np.isfinite(xvalue)]),np.max(xvalue[np.isfinite(xvalue)]),0.1)
    ax.set_xscale('log')
    ax.plot(xaxis, 0*xaxis,'k')
    medoffset = np.nanmedian(residual[np.isfinite(residual)])
    ax.plot(xaxis, 0*xaxis+medoffset,'k-.')
    
#manager = plt.get_current_fig_manager()
#manager.window.showMaximized()
fig.tight_layout()

fig, axes = plt.subplots(2,2, figsize = (13.5,6.5))
fig.suptitle("NSA vs. Reines")
for i in range(len(numlines)):
    
    numline = numlines[i]
    denline = denlines[i]
    ax = axes.ravel()[numlines.index(numline)]
    if numline == 'sii_6731_flux':        
        reinesnum = (reines[reineslines[line]].loc[nsa['NSAID'][ndx]]+\
                          reines['SII-6716'].loc[nsa['NSAID'][ndx]]) 
    else:
        reinesnum = reines[reineslines[numline]].loc[nsa['NSAID'][ndx]]
    reinesden = reines[reineslines[denline]].loc[nsa['NSAID'][ndx]]
    xvalue = 10**np.log10(reinesnum/reinesden)
    print(numline)
    nsanum = nsa[nsalines[numline]][ndx]
    nsaden = nsa[nsalines[denline]][ndx]
    yvalue = 10**np.log10(nsanum/nsaden)
    residual = (yvalue - xvalue)/xvalue
    ax.plot(xvalue, residual, 'r.')
    ax.set_xlabel('Reines '+numline+'/'+denline)
    ax.set_ylabel('(NSA - Reines)/Reines')#+numline+'/'+denline)
    ax.set_ylim(-0.2,1)
    ax.set_xscale('log')
    xaxis = np.arange(np.min(xvalue[np.isfinite(xvalue)]),np.max(xvalue[np.isfinite(xvalue)]),0.1)
    ax.plot(xaxis, 0*xaxis,'k')
    medoffset = np.nanmedian(residual[np.isfinite(residual)])
    ax.plot(xaxis, 0*xaxis+medoffset,'k-.')
    print(np.nanmedian(residual[np.isfinite(residual)]))
    print(min(residual), max(residual))
#manager = plt.get_current_fig_manager()
#manager.window.showMaximized()
fig.tight_layout()
