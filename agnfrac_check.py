# -*- coding: utf-8 -*-
"""
Created on Tue May 12 12:21:06 2020

@author: mugdhapolimera

Checking AGN fraction as a function of environment in different Declination 
bins
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.stats import binom_conf_interval
from scipy.io.idl import readsav
resfull = pd.read_csv('RESOLVE_snr5_master_new.csv')
#resflag = pd.read_csv('resolve_emlineclass_full_snr5_master_new.csv')
#resconf = pd.read_csv('RESOLVE_snr5_master_conf_bary.csv')
resfullbary = pd.read_csv('RESOLVE_full_dext_snr5_jhu.csv')
resexcluded = [x for x in list(resfull.name) if x not in list(resfullbary.name)]

resfull = resfullbary
resflag = pd.read_csv('resolve_emlineclass_dext_snr5_jhu.csv')
resflagagn = resflag.sftoagn | resflag.agntosf | resflag.defagn | resflag.composite
resconf = pd.DataFrame({'name': resflag.galname,
                        'confidence_level': resflagagn*1.0,
                        'jhu' : resflagagn*1.0,
                        'port': np.zeros(len(resflag.galname))})
resconf.confidence_level[resconf.confidence_level == 0] = -1
#resconf = pd.read_csv('RESOLVE_snr5_master_conf_bary.csv')

resfull['logmbary'] = np.log10(10**resfull.logmstar + 10**resfull.logmgas)
resfull.index = resfull.name
resflag.index = resflag.galname
resconf.index = resconf.name
#resconf.confidence_level = resconf.jhu + resconf.port

rescatphot = readsav('resolvecatalogphot.dat')
rescat = readsav('resolvecatalog.dat')
resdat_ndx = [x for x in range(len(rescat.name)) if list(rescat.name)[x] in \
           list(resfull.name) ]
resfull['mur'] = np.zeros(len(resfull.name))
resfull.loc[list(rescat.name[resdat_ndx]),'mur']= np.array(rescatphot.mur_radr23p75[resdat_ndx]).byteswap().newbyteorder()
#resflag['defstarform'][resconf.confidence_level < 1] = True
#resflag['sftoagn'][resconf.confidence_level < 1] = False
#resflag['agntosf'][resconf.confidence_level < 1] = False
#resflag['composite'][resconf.confidence_level < 1] = False
#resflag['defagn'][resconf.confidence_level < 1] = False
#
#res = resfull
resndx = (resfull.logmbary > 9.2) # | (resfull.absrmag < -17.33) 
res = resfull[resndx]
resflag = resflag[resndx]
resconf = resconf.loc[res.name]
ecofull = pd.read_csv('ECO_snr5_master_new.csv')
#ecofullflag = pd.read_csv('eco_emlineclass_full_snr5_master_new.csv')
#ecoconf = pd.read_csv('ECO_snr5_master_conf_bary.csv')
ecofullbary = pd.read_csv('ECO_full_dext_snr5_jhu.csv')
ecoexcluded = [x for x in list(ecofull.name) if x not in list(ecofullbary.name)]

ecofull = ecofullbary
ecofullflag = pd.read_csv('eco_emlineclass_dext_snr5_jhu.csv')
#ecoconf = pd.read_csv('ECO_snr5_master_conf_bary.csv')
ecofull.index = ecofull.name
ecofullflag.index = ecofullflag.galname
ecofull['logmbary'] = np.log10(10**ecofull.logmstar + 10**ecofull.logmgas)
ecocat = readsav('eco_wresa_032918.dat')
ecodat_ndx = [x for x in range(len(ecocat.econames)) if list(ecocat.econames)[x] in \
           list(ecofull.name) ]
ecofull['mur'] = np.zeros(len(ecofull.name))
ecofull.loc[list(ecocat.econames[ecodat_ndx]),'mur']= np.array(ecocat['mur_radr23p75'][ecodat_ndx]).byteswap().newbyteorder()


econdx = (((ecofull.cz < 4500) | (ecofull.dedeg < 0) | (ecofull.dedeg > 5))) #| \
#         (ecofull.resname == 'notinresolve') )
#eco = ecofull[econdx]
#ecoflag = ecofullflag[econdx]
#ecoconf.index = ecoconf.name
#ecoconf.confidence_level = ecoconf.jhu + ecoconf.port
eco = ecofull[econdx & ((ecofull.logmbary > 9.2))]# | (ecofull.absrmag < -17.33))] 
ecoflag = ecofullflag[econdx & ((ecofull.logmbary > 9.2))]# | (ecofull.absrmag < -17.33))]
#ecoconf = ecoconf.loc[eco.name]
ecoflagagn = ecoflag.sftoagn | ecoflag.agntosf | ecoflag.defagn | ecoflag.composite
ecoconf = pd.DataFrame({'name': ecoflag.galname,
                        'confidence_level': ecoflagagn*1.0,
                        'jhu' : ecoflagagn*1.0,
                        'port': np.zeros(len(ecoflag.galname))})
ecoconf.confidence_level[ecoconf.confidence_level == 0] = -1

#ecoflag['defstarform'][ecoconf.confidence_level < 1] = True
#ecoflag['sftoagn'][ecoconf.confidence_level < 1] = False
#ecoflag['agntosf'][ecoconf.confidence_level < 1] = False
#ecoflag['composite'][ecoconf.confidence_level < 1] = False
#ecoflag['defagn'][ecoconf.confidence_level < 1] = False
#
ecog3full = pd.read_csv('ecog3.csv')
ecog3full.index = ecog3full.name
ecog3ndx = (((ecog3full.cz < 4500) | (ecog3full.dedeg < 0) | (ecog3full.dedeg > 5)) | \
         (ecog3full.resname == 'notinresolve'))
ecog3 = ecog3full[ecog3ndx]

#eco_lssdens_dat = readsav('eco_lss_dens_mhalo_N=3_dv500_091016.dat').lss_dens

#eco.logmh = ecog3.loc[ecoflag.galname].logMh_g3
#eco.grp = ecog3.loc[ecoflag.galname].grp_g3
#eco.grpn = ecog3.loc[ecoflag.galname].grpn_g3
#eco.fc = ecog3.loc[ecoflag.galname].central_g3


def agnfrac_env(galtype = 'all'):
    if galtype == 'all':
        res_sub_ndx = res.logmstar > 0
        eco_sub_ndx = eco.logmstar > 0
    elif galtype == 'dwarf':
        res_sub_ndx = (res.logmstar < 9.5)
        eco_sub_ndx = (eco.logmstar < 9.5)
    
    res_sub = res[res_sub_ndx]
    res_sub_flag = resflag[res_sub_ndx]
    res_sub_conf = resconf[res_sub_ndx]
    eco_sub = eco[eco_sub_ndx]
    eco_sub_flag = ecoflag[eco_sub_ndx]
    eco_sub_conf = ecoconf[eco_sub_ndx]
    
    res_sub_agn = res_sub_flag.sftoagn |  res_sub_flag.defagn | \
                    res_sub_flag.composite | res_sub_flag.agntosf
    eco_sub_agn = eco_sub_flag.sftoagn |  eco_sub_flag.defagn | \
                    eco_sub_flag.composite | eco_sub_flag.agntosf
    res_sub_agn = res_sub_conf.confidence_level >= 0
    eco_sub_agn = eco_sub_conf.confidence_level >= 0
    #res_sub_agnfrac = 100.0*np.sum(res_sub_agn)/len(res_sub)
    #eco_sub_agnfrac = 100.0*np.sum(eco_sub_agn)/len(eco_sub)

    
    bins = np.arange(10.5, 15.5, 0.5)
    res_sub_env = np.zeros(len(bins) - 1)
    eco_sub_env = np.zeros(len(bins) - 1)
    res_sub_agn_env = np.zeros(len(bins) - 1)
    eco_sub_agn_env = np.zeros(len(bins) - 1)
    mh_bins = np.zeros(len(bins) - 1)
#    for i in range(10):
#        lowmh = 10.5 + i*1
#        highmh = lowmh + 0.5
#        
#        mh_bins[i] = highmh
#        res_env_ndx = (res_sub.logmh > lowmh) & (res_sub.logmh < highmh)
#        res_sub_env[i] = len(res_sub[res_env_ndx])
#        res_sub_agn_env[i] = np.sum(res_sub_agn & res_env_ndx)
#        
#        eco_env_ndx = (eco_sub.logmh > lowmh) & (eco_sub.logmh < highmh)
#        eco_sub_env[i] = len(eco_sub[eco_env_ndx])
#        eco_sub_agn_env[i] = np.sum(eco_sub_agn & eco_env_ndx)
#    
    
    res_sub_env = np.histogram(res_sub.logmh, bins = bins)[0]
    res_sub_agn_env = np.histogram(res_sub[res_sub_agn].logmh, bins = bins)[0]
    eco_sub_env = np.histogram(eco_sub.logmh, bins = bins)[0]
    eco_sub_agn_env = np.histogram(eco_sub[eco_sub_agn].logmh, bins = bins)[0]
    
    res_sub_agnfrac = 100.0*(res_sub_agn_env)/(res_sub_env)
    eco_sub_agnfrac = 100.0*(eco_sub_agn_env)/(eco_sub_env)

    reslow = np.zeros(len(res_sub_agnfrac))
    resup = np.zeros(len(res_sub_agnfrac))
    ecolow = np.zeros(len(res_sub_agnfrac))
    ecoup = np.zeros(len(res_sub_agnfrac))
    for i in range(len(res_sub_agn_env)):
        if res_sub_env[i] > 0:
        
            reslowlim, resuplim = 100*binom_conf_interval(res_sub_agn_env[i],\
                                                      res_sub_env[i])
            resup[i]= resuplim - res_sub_agnfrac[i]
            reslow[i] = res_sub_agnfrac[i] - reslowlim 
        else:
#            res_sub_agnfrac[i] = 0
            resup[i]= 0
            reslow[i] = 0
    for i in range(len(res_sub_agn_env)):
        if eco_sub_env[i] > 0:

            ecolowlim, ecouplim = 100*binom_conf_interval(eco_sub_agn_env[i],\
                                                          eco_sub_env[i])
            ecoup[i]= ecouplim - eco_sub_agnfrac[i]
            ecolow[i] = eco_sub_agnfrac[i] - ecolowlim 
        else:
#            eco_sub_agnfrac[i] = 0
            ecoup[i]= 0
            ecolow[i] = 0
    
#    res_sub_agnfrac = res_sub_agn_env
#    eco_sub_agnfrac = eco_sub_agn_env
#    ecoup= ecoup/100*len(eco_sub)
#    ecolow = ecolow/100*len(eco_sub)
#    resup= ecoup/100*len(res_sub)
#    reslow = ecolow/100*len(res_sub)
#    
    
    return res_sub_agnfrac, reslow, resup, \
            eco_sub_agnfrac, ecolow, ecoup,mh_bins 

def agnfrac_dens(galtype = 'all', columnname = 'den1mpc'):
    
    if columnname == 'den1mpc':
        if galtype == 'all':
            eco_sub_ndx = (eco.den1mpc > -99) & (eco.logmstar > 0)
            res_sub_ndx = ~econdx & (ecofull.logmstar > 0) & (ecofull.den1mpc > -99)
        elif galtype == 'dwarf':
            res_sub_ndx = ~econdx & (ecofull.logmstar < 9.5) & (ecofull.den1mpc > -99)
            eco_sub_ndx = (eco.logmstar < 9.5) & (eco.den1mpc > -99)
        
    else: #if columnname == 'dens_s':
        if galtype == 'all':
            eco_sub_ndx = (eco.dens_s > 0) & (eco.logmstar > 0)
            res_sub_ndx = ~econdx & (ecofull.logmstar > 0) & (ecofull.dens_s > 0)
        elif galtype == 'dwarf':
            res_sub_ndx = (res.logmstar < 9.5) & (res[columnname] > 0)
            eco_sub_ndx = (eco.logmstar < 9.5) & (eco[columnname] > 0)
    
    res_sub = res[res_sub_ndx]
    res_sub_flag = resflag[res_sub_ndx]
    res_sub_conf = resconf[res_sub_ndx]
    eco_sub = eco[eco_sub_ndx]
    eco_sub_flag = ecoflag[eco_sub_ndx]
    eco_sub_conf = ecoconf[eco_sub_ndx]
    
    res_sub_agn = res_sub_flag.sftoagn |  res_sub_flag.defagn | \
                    res_sub_flag.composite | res_sub_flag.agntosf
    eco_sub_agn = eco_sub_flag.sftoagn |  eco_sub_flag.defagn | \
                    eco_sub_flag.composite | eco_sub_flag.agntosf
    res_sub_agn = res_sub_conf.confidence_level >= 0
    eco_sub_agn = eco_sub_conf.confidence_level >= 0

    bins = np.arange(int(min(res_sub[columnname])), \
                     int(max(res_sub[columnname])), 0.2)
#    res_sub_env = np.zeros(len(bins) - 1)
#    eco_sub_env = np.zeros(len(bins) - 1)
#    res_sub_agn_env = np.zeros(len(bins) - 1)
#    eco_sub_agn_env = np.zeros(len(bins) - 1)
    mh_bins = bins[:-1]#np.zeros(len(bins) - 1)
#    for i in range(10):
#        lowmh = 10.5 + i*1
#        highmh = lowmh + 0.5
#        
#        mh_bins[i] = highmh
#        res_env_ndx = (res_sub.logmh > lowmh) & (res_sub.logmh < highmh)
#        res_sub_env[i] = len(res_sub[res_env_ndx])
#        res_sub_agn_env[i] = np.sum(res_sub_agn & res_env_ndx)
#        
#        eco_env_ndx = (eco_sub.logmh > lowmh) & (eco_sub.logmh < highmh)
#        eco_sub_env[i] = len(eco_sub[eco_env_ndx])
#        eco_sub_agn_env[i] = np.sum(eco_sub_agn & eco_env_ndx)
#    
    
    res_sub_env = np.histogram(res_sub[columnname], bins = bins)[0]
    res_sub_agn_env = np.histogram(res_sub[res_sub_agn][columnname], bins = bins)[0]
    eco_sub_env = np.histogram(eco_sub[columnname], bins = bins)[0]
    eco_sub_agn_env = np.histogram(eco_sub[eco_sub_agn][columnname], bins = bins)[0]
    
    print(res_sub_env, res_sub_agn_env)
    print(eco_sub_env, eco_sub_agn_env)
    res_sub_agnfrac = 100.0*(res_sub_agn_env)/(res_sub_env)
    eco_sub_agnfrac = 100.0*(eco_sub_agn_env)/(eco_sub_env)

    reslow = np.zeros(len(res_sub_agnfrac))
    resup = np.zeros(len(res_sub_agnfrac))
    ecolow = np.zeros(len(res_sub_agnfrac))
    ecoup = np.zeros(len(res_sub_agnfrac))
    for i in range(len(res_sub_agn_env)):
        if res_sub_env[i] > 0:
        
            reslowlim, resuplim = 100*binom_conf_interval(res_sub_agn_env[i],\
                                                      res_sub_env[i])
            resup[i]= resuplim - res_sub_agnfrac[i]
            reslow[i] = res_sub_agnfrac[i] - reslowlim 
        else:
#            res_sub_agnfrac[i] = 0
            resup[i]= 0
            reslow[i] = 0
    for i in range(len(res_sub_agn_env)):
        if eco_sub_env[i] > 0:

            ecolowlim, ecouplim = 100*binom_conf_interval(eco_sub_agn_env[i],\
                                                          eco_sub_env[i])
            ecoup[i]= ecouplim - eco_sub_agnfrac[i]
            ecolow[i] = eco_sub_agnfrac[i] - ecolowlim 
        else:
#            eco_sub_agnfrac[i] = 0
            ecoup[i]= 0
            ecolow[i] = 0
    
#    res_sub_agnfrac = res_sub_agn_env
#    eco_sub_agnfrac = eco_sub_agn_env
#    ecoup= ecoup/100*len(eco_sub)
#    ecolow = ecolow/100*len(eco_sub)
#    resup= ecoup/100*len(res_sub)
#    reslow = ecolow/100*len(res_sub)
#    
    
    return res_sub_agnfrac, reslow, resup, \
            eco_sub_agnfrac, ecolow, ecoup,mh_bins 


ecoparent = pd.read_csv('ECO_inobssample.csv')
ecoparent.index = ecoparent.name
ecoparentndx = ((ecoparent.dedeg < 0) | (ecoparent.dedeg > 5) & \
         (ecoparent.resname == 'notinresolve'))
ecoparent = ecoparent[ecoparentndx]
resparent = pd.read_csv('RESOLVE_inobssample.csv')
resparent.index = resparent.name
bins = np.arange(10.5,15.5,0.5)
fig = plt.figure()
ax1 = plt.subplot(121)
ax1.hist(ecoparent.logmh, bins = bins, histtype = 'step', 
         edgecolor = 'k', lw = 3, label = 'ECO Parent Sample (excluding RESOLVE-A)')#,normed = True)
ax1.hist(eco.logmh, bins = bins, histtype = 'step', hatch = 'x',
         edgecolor = 'orange', lw = 2, label = 'ECO SELs (excluding RESOLVE-A)')#, normed = True)
#ax1.hist(res.logmh, bins = bins, histtype = 'step', hatch = 'x',
#         edgecolor = 'orange', lw = 2, ls = '--',
#         label = 'RESOLVE SELs')#, normed = True)
#ax1.hist(resparent.logmh, bins = bins, histtype = 'step',  
#         edgecolor = 'black', lw = 3, ls = '--',
#         label = 'RESOLVE Parent Sample')#, normed = True)
ax1.legend()
ax1.set_xlabel('log(M$_{halo}$/M$_{\odot}$)', fontsize = 15) 
ax1.set_ylabel('Number of galaxies', fontsize = 15)
ax2 = plt.subplot(122)
ax2.hist(resparent.logmh, bins = bins, histtype = 'step',  
         edgecolor = 'black', lw = 3, label = 'RESOLVE Parent Sample')#, normed = True)
ax2.hist(res.logmh, bins = bins, histtype = 'step', hatch = 'x',
         edgecolor = 'orange', lw = 2, label = 'RESOLVE SELs')#, normed = True)
ax2.legend()
ax2.set_xlabel('log(M$_{halo}$/M$_{\odot}$)', fontsize = 15)
ax2.set_ylabel('Number of galaxies', fontsize = 15)


ecodwarf = eco[eco.logmstar < 9.5]
resdwarf = res[res.logmstar < 9.5]

ecodwarf_hist = np.histogram(ecodwarf.logmh, bins = bins)
resdwarf_hist = np.histogram(resdwarf.logmh, bins = bins)
#res_agn = resflag.sftoagn |  resflag.defagn | \
#                    resflag.composite | resflag.agntosf
#eco_agn = ecoflag.sftoagn |  ecoflag.defagn | \
#                    ecoflag.composite | ecoflag.agntosf
res_agn = resconf.confidence_level >= 0
eco_agn = ecoconf.confidence_level >= 0
ecogiantagn = eco[(eco.logmstar > 9.5) & (eco_agn)]
resgiantagn = res[(res.logmstar > 9.5) & (res_agn)]
ecodwarfagn = ecodwarf[eco_agn]
resdwarfagn = resdwarf[res_agn]
ecoagn = eco[eco_agn]
resagn = res[res_agn]
resexcludedagn = [x for x in list(resagn.name) \
                  if (x not in list(resfullbary.name)) ]
ecoexcludedagn = [x for x in list(ecoagn.name) \
                  if (x not in list(ecofullbary.name)) & (x not in list(eco.name)) ]

resjhuandport = (resconf.confidence_level == 2) | \
            (resconf.confidence_level == -2) |(resconf.confidence_level == 0)
resjhuandportname_agn = list(resconf.loc[resdwarfagn.name][resconf.confidence_level == 2].name)
resjhuandportname = list(resconf.loc[resdwarf.name][resjhuandport].name)

ecojhuandport = (ecoconf.confidence_level == 2) | \
            (ecoconf.confidence_level == -2) |(ecoconf.confidence_level == 0)
ecojhuandportname_agn = list(ecoconf.loc[ecodwarfagn.name][ecoconf.confidence_level == 2].name)
ecojhuandportname = list(ecoconf.loc[ecodwarf.name][ecojhuandport].name)

#resdwarf = resdwarf.loc[resjhuandportname]
#resdwarfagn = resdwarfagn.loc[resjhuandportname_agn]
#
#ecodwarf = ecodwarf.loc[ecojhuandportname]
#ecodwarfagn = ecodwarfagn.loc[ecojhuandportname_agn]

ecodwarfagn_hist = np.histogram(ecodwarfagn.logmh, bins = bins)
resdwarfagn_hist = np.histogram(resdwarfagn.logmh, bins = bins)

econonisodwarf = ecodwarf[ecodwarf.grpn > 1]
ecoisodwarf = ecodwarf[ecodwarf.grpn == 1]
ecoisodwarfagn = ecodwarfagn[ecodwarfagn.grpn == 1]
econonisodwarfagn = ecodwarfagn[ecodwarfagn.grpn > 1]
resnonisodwarfagn = resdwarfagn[resdwarfagn.grpn > 1]
resisodwarfagn = resdwarfagn[resdwarfagn.grpn == 1]
resisodwarf = resdwarf[resdwarf.grpn == 1]
resnonisodwarf = resdwarf[resdwarf.grpn > 1]
ecocendwarfagn = ecodwarfagn[ecodwarfagn.fc == 1]
rescendwarfagn = resdwarfagn[resdwarfagn.fc == 1]
ressatdwarfagn = resdwarfagn[resdwarfagn.fc == 0]
ressatdwarf = resdwarf[resdwarf.fc == 0]
rescendwarf = resdwarf[resdwarf.fc == 1]
ecosatdwarfagn = ecodwarfagn[ecodwarfagn.fc == 0]
ecosatdwarf = ecodwarf[ecodwarf.fc == 0]
ecocendwarf = ecodwarf[ecodwarf.fc == 1]

##ECO Dwarf AGN % breakdown by type
ecodwarfagnpc = len(ecodwarfagn)*100.0/len(ecodwarf)
ecodwarfagnpc_low, ecodwarfagnpc_up = 100*binom_conf_interval(len(ecodwarfagn),len(ecodwarf)) - ecodwarfagnpc

ecoisodwarfagnpc = len(ecoisodwarfagn)*100.0/len(ecoisodwarf)
ecoisodwarfagnpc_low, ecoisodwarfagnpc_up = 100*binom_conf_interval(len(ecoisodwarfagn),
                                            len(ecoisodwarf)) - ecoisodwarfagnpc

econonisodwarfagnpc = len(econonisodwarfagn)*100.0/len(econonisodwarf)
econonisodwarfagnpc_low, econonisodwarfagnpc_up = 100*binom_conf_interval(len(econonisodwarfagn),
                                            len(econonisodwarf)) - econonisodwarfagnpc

ecocendwarfagnpc = len(ecocendwarfagn)*100.0/len(ecocendwarf)
ecocendwarfagnpc_low, ecocendwarfagnpc_up = 100*binom_conf_interval(len(ecocendwarfagn),
                                          len(ecocendwarf)) - ecocendwarfagnpc

ecosatdwarfagnpc = len(ecosatdwarfagn)*100.0/len(ecosatdwarf)
ecosatdwarfagnpc_low, ecosatdwarfagnpc_up = 100*binom_conf_interval(len(ecosatdwarfagn),
                                          len(ecosatdwarf)) - ecosatdwarfagnpc

#                                                                    
#ecodwarfagnpc = np.sum(ecodwarfagn.ccr)*100.0/np.sum(ecodwarf.ccr)
#ecodwarfagnpc_low, ecodwarfagnpc_up = 100*binom_conf_interval(\
#                        np.sum(ecodwarfagn.ccr),np.sum(ecodwarf.ccr)) - ecodwarfagnpc
#
#ecoisodwarfagnpc = np.sum(ecoisodwarfagn.ccr)*100.0/np.sum(ecoisodwarf.ccr)
#ecoisodwarfagnpc_low, ecoisodwarfagnpc_up = 100*binom_conf_interval(\
#                        np.sum(ecoisodwarfagn.ccr),np.sum(ecoisodwarf.ccr)) - ecoisodwarfagnpc
#
#econonisodwarfagnpc = np.sum(econonisodwarfagn.ccr)*100.0/np.sum(econonisodwarf.ccr)
#econonisodwarfagnpc_low, econonisodwarfagnpc_up = 100*binom_conf_interval(\
#                        np.sum(econonisodwarfagn.ccr),np.sum(econonisodwarf.ccr)) - econonisodwarfagnpc
#
#ecocendwarfagnpc = np.sum(ecocendwarfagn.ccr)*100.0/np.sum(ecocendwarf.ccr)
#ecocendwarfagnpc_low, ecocendwarfagnpc_up = 100*binom_conf_interval(\
#                        np.sum(ecocendwarfagn.ccr),np.sum(ecocendwarf.ccr)) - ecocendwarfagnpc
#
#ecosatdwarfagnpc = np.sum(ecosatdwarfagn.ccr)*100.0/np.sum(ecosatdwarf.ccr)
#ecosatdwarfagnpc_low, ecosatdwarfagnpc_up = 100*binom_conf_interval(\
#                        np.sum(ecosatdwarfagn.ccr),np.sum(ecosatdwarf.ccr)) - ecosatdwarfagnpc
#                                                                    
#                                                                    
##RESOLVE Dwarf AGN % breakdown by type
resdwarfagnpc = len(resdwarfagn)*100.0/len(resdwarf)
resdwarfagnpc_low, resdwarfagnpc_up = 100*binom_conf_interval(len(resdwarfagn),len(resdwarf)) - resdwarfagnpc

resisodwarfagnpc = len(resisodwarfagn)*100.0/len(resisodwarf)
resisodwarfagnpc_low, resisodwarfagnpc_up = 100*binom_conf_interval(len(resisodwarfagn),
                                            len(resisodwarf)) - resisodwarfagnpc

resnonisodwarfagnpc = len(resnonisodwarfagn)*100.0/len(resnonisodwarf)
resnonisodwarfagnpc_low, resnonisodwarfagnpc_up = 100*binom_conf_interval(len(resnonisodwarfagn),
                                            len(resnonisodwarf)) - resnonisodwarfagnpc

rescendwarfagnpc = len(rescendwarfagn)*100.0/len(rescendwarf)
rescendwarfagnpc_low, rescendwarfagnpc_up = 100*binom_conf_interval(len(rescendwarfagn),
                                          len(rescendwarf)) - rescendwarfagnpc

ressatdwarfagnpc = len(ressatdwarfagn)*100.0/len(ressatdwarf)
ressatdwarfagnpc_low, ressatdwarfagnpc_up = 100*binom_conf_interval(len(ressatdwarfagn),
                                          len(ressatdwarf)) - ressatdwarfagnpc
#
#resdwarfagnpc = np.sum(resdwarfagn.ccr)*100.0/np.sum(resdwarf.ccr)
#resdwarfagnpc_low, resdwarfagnpc_up = 100*binom_conf_interval(\
#                        np.sum(resdwarfagn.ccr),np.sum(resdwarf.ccr)) - resdwarfagnpc
#                                                                    
#resisodwarfagnpc = np.sum(resisodwarfagn.ccr)*100.0/np.sum(resisodwarf.ccr)
#resisodwarfagnpc_low, resisodwarfagnpc_up = 100*binom_conf_interval(\
#                        np.sum(resisodwarfagn.ccr),np.sum(resisodwarf.ccr)) - resisodwarfagnpc
#
#resnonisodwarfagnpc = np.sum(resnonisodwarfagn.ccr)*100.0/np.sum(resnonisodwarf.ccr)
#resnonisodwarfagnpc_low, resnonisodwarfagnpc_up = 100*binom_conf_interval(\
#                        np.sum(resnonisodwarfagn.ccr),np.sum(resnonisodwarf.ccr)) - resnonisodwarfagnpc
#
#rescendwarfagnpc = np.sum(rescendwarfagn.ccr)*100.0/np.sum(rescendwarf.ccr)
#rescendwarfagnpc_low, rescendwarfagnpc_up = 100*binom_conf_interval(\
#                        np.sum(rescendwarfagn.ccr),np.sum(rescendwarf.ccr)) - rescendwarfagnpc
#
#ressatdwarfagnpc = np.sum(ressatdwarfagn.ccr)*100.0/np.sum(ressatdwarf.ccr)
#ressatdwarfagnpc_low, ressatdwarfagnpc_up = 100*binom_conf_interval(\
#                        np.sum(ressatdwarfagn.ccr),np.sum(ressatdwarf.ccr)) - ressatdwarfagnpc
                                                                    
xaxis = np.arange(1,11,2)
respc = [resdwarfagnpc, resisodwarfagnpc, resnonisodwarfagnpc, rescendwarfagnpc, \
         ressatdwarfagnpc]
ecopc = [ecodwarfagnpc, ecoisodwarfagnpc, econonisodwarfagnpc, ecocendwarfagnpc, \
         ecosatdwarfagnpc]

respc_low = [resdwarfagnpc_low, resisodwarfagnpc_low, resnonisodwarfagnpc_low, rescendwarfagnpc_low, \
         ressatdwarfagnpc_low]
ecopc_low = [ecodwarfagnpc_low, ecoisodwarfagnpc_low, econonisodwarfagnpc_low, ecocendwarfagnpc_low, \
         ecosatdwarfagnpc_low]

respc_up = [resdwarfagnpc_up, resisodwarfagnpc_up, resnonisodwarfagnpc_up, rescendwarfagnpc_up, \
         ressatdwarfagnpc_up]
ecopc_up = [ecodwarfagnpc_up, ecoisodwarfagnpc_up, econonisodwarfagnpc_up, ecocendwarfagnpc_up, \
         ecosatdwarfagnpc_up]
                                                                    
plt.figure()
#plt.suptitle('Mass limited sample; AGN w/ conf = 2')
plt.errorbar(xaxis, respc, fmt= 'bs', yerr = [-1.0*np.array(respc_low), respc_up], label = 'RESOLVE')
plt.xticks(xaxis, ['Overall', 'Iso. Dwarf', 'Non-iso. Dwarf', 'Cen. Dwarf', 'Sat. Dwarf'])
plt.errorbar(xaxis+0.2, ecopc, fmt= 'rs', yerr = [-1.0*np.array(ecopc_low), ecopc_up], label = 'ECO (excluding RESOLVE-A)')
plt.ylabel('AGN % in Dwarf SELs')
plt.legend()
                                                                    
plt.figure()
plt.hist(econonisodwarf.logmh, bins = bins, label = 'ECO non-iso dwarf SELs',
         histtype = 'step', normed = True)
plt.hist(resnonisodwarf.logmh, bins = bins, label = 'RESOLVE non-iso dwarf SELs',
         histtype = 'step', normed = True)
#plt.hist(ecoisodwarf.logmh, bins = bins, label = 'ECO iso dwarf SELs',
#         histtype = 'step', normed = True)
#plt.hist(resisodwarf.logmh, bins = bins, label = 'RESOLVE iso dwarf SELs',
#         histtype = 'step', normed = True)
#plt.hist(ecoisodwarfagn.logmh, bins = bins, label = 'ECO iso dwarf SEL AGN',
#         histtype = 'step', normed = True)
#plt.hist(resisodwarfagn.logmh, bins = bins, label = 'RESOLVE iso dwarf SEL AGN',
#         histtype = 'step', normed = True)
plt.hist(econonisodwarfagn.logmh, bins = bins, label = 'ECO non-iso dwarf SEL AGN',
         histtype = 'step', normed = True)
plt.hist(resnonisodwarfagn.logmh, bins = bins, label = 'RESOLVE non-iso dwarf SEL AGN',
         histtype = 'step', normed = True)

plt.legend()

plt.figure()
#c = ['g','b',r','m','k']
res_sub_agnfrac, reslow, resup, \
        eco_sub_agnfrac, ecolow, ecoup,mhbins = \
        agnfrac_env('dwarf')
reslow[reslow<0.1] = 0
ecolow[ecolow<0.1] = 0

plt.bar(bins[:9]+0.25,eco_sub_agnfrac,yerr = [ecolow,ecoup],color = 'none',
        edgecolor = 'blue', linewidth = 3, label = 'ECO SEL Dwarfs (excluding RESOLVE-A)',
        align = 'center',width = np.diff(bins)[0],capsize = 5, hatch = 'x',
        ecolor = 'blue')        
plt.bar(bins[:9]+0.25,res_sub_agnfrac,yerr = [reslow,resup],color = 'none',
        edgecolor = 'orange', linewidth = 3, label = 'RESOLVE SEL Dwarfs',
        ecolor = 'orange',
        align = 'center',width = np.diff(bins)[0],capsize = 5)        
#plt.errorbar(bins[:9]+0.25,eco_sub_agnfrac,fmt = '.',yerr = [ecolow,ecoup],\
#             label = 'ECO SEL Dwarfs')
#plt.errorbar(bins[:9]+0.25,res_sub_agnfrac,fmt = '.',yerr = [reslow,resup],\
#             label = 'RESOLVE SEL Dwarfs')

#plt.errorbar(bins[:9]+0.5,eco_sub_agnfrac,marker = '.', drawstyle = 'steps-mid', 
#             yerr = [ecolow,ecoup],label = 'ECO SEL Dwarfs')
#plt.errorbar(bins[:9]-0.1+0.5,res_sub_agnfrac,marker = '.', drawstyle = 'steps-mid', 
#             yerr = [reslow,resup],label = 'RESOLVE SEL Dwarfs')

plt.legend(loc = 'upper left')
plt.xlabel('Log(M$_{halo}$/M$_\odot$)', fontsize = 15)
plt.ylabel('Dwarf AGN Frequency', fontsize = 15)
#plt.ylim(0.,60)
plt.xlim(10.25,15.25)

plt.figure()
res_sub_agnfrac, reslow, resup, \
        eco_sub_agnfrac, ecolow, ecoup,mhbins = \
        agnfrac_env('dwarf')
reslow[reslow<0.1] = 0
ecolow[ecolow<0.1] = 0

plt.errorbar(bins[:9]+0.25,eco_sub_agnfrac,fmt = '.-',yerr = [ecolow,ecoup],\
             label = 'ECO SEL Dwarfs (excluding RESOLVE-A)', capsize = 5)
plt.errorbar(bins[:9]+0.25,res_sub_agnfrac,fmt = '.-',yerr = [reslow,resup],\
             label = 'RESOLVE SEL Dwarfs', capsize = 5)

#plt.errorbar(bins[:9]+0.25,eco_sub_agnfrac,marker = '.', drawstyle = 'steps-mid', 
#             yerr = [ecolow,ecoup],label = 'ECO SEL Dwarfs')
#plt.errorbar(bins[:9]+0.25,res_sub_agnfrac,marker = '.', drawstyle = 'steps-mid', 
#             yerr = [reslow,resup],label = 'RESOLVE SEL Dwarfs')

plt.legend(loc = 'upper left')
plt.xlabel('Log(M$_{halo}$/M$_\odot$)', fontsize = 15)
plt.ylabel('Dwarf AGN Frequency', fontsize = 15)
#plt.ylim(0.,60)
plt.xlim(10.5,15.0)
plt.xticks(np.arange(10.5,15.5,0.5))

#Den1MPC
#fig = plt.figure()
#ecoden = eco.den1mpc[(eco.logmstar < 9.5) & (eco.den1mpc > -99)]
#ecomh = eco.logmh[(eco.logmstar < 9.5) & (eco.den1mpc > -99)]
#resaden = ecofull.den1mpc[~econdx & (ecofull.logmstar < 9.5) & (ecofull.den1mpc > -99)]
#resamh = ecofull.logmh[~econdx & (ecofull.logmstar < 9.5) & (ecofull.den1mpc > -99)]
#ecodwarfndx = (eco.logmstar < 9.5) & (eco.den1mpc > -99)
##mstarxaxis, mstarkde = ppty_kde(ecoden, xlabel = 'normalized environmental density smoothed at 1Mpc scale', plttype = ptype,
##                                label = 'ECO SEL')
##mstarxaxis, mstarkde = ppty_kde(resaden, xlabel = 'normalized environmental density smoothed at 1Mpc scale',plttype = ptype,
##                                label = 'RESOLVE SEL')
##plt.legend()
#
##kde_errorband(np.array(ecoden),-1,5,'red', 'ECO SEL Dwarfs (excluding RESOLVE-A)', \
##              'normalized environmental density smoothed at 1Mpc scale')
##kde_errorband(np.array(resaden),-1,5,'lime', 'RESOLVE-A SEL Dwarfs', \
##              'normalized environmental density smoothed at 1Mpc scale')
#
plt.figure()
#c = ['g','b',r','m','k']
res_sub_agnfrac, reslow, resup, \
        eco_sub_agnfrac, ecolow, ecoup,mhbins = \
        agnfrac_dens('dwarf', 'logmbary')
reslow[reslow<0.1] = 0
ecolow[ecolow<0.1] = 0

#plt.bar(mhbins+0.25,eco_sub_agnfrac,yerr = [ecolow,ecoup],color = 'none',
#        edgecolor = 'blue', linewidth = 3, label = 'ECO SEL Dwarfs (excluding RESOLVE-A)',
#        align = 'center',width = np.diff(bins)[0],capsize = 5, hatch = 'x',
#        ecolor = 'blue')        
#plt.bar(mhbins+0.25,res_sub_agnfrac,yerr = [reslow,resup],color = 'none',
#        edgecolor = 'orange', linewidth = 3, label = 'RESOLVE SEL Dwarfs',
#        ecolor = 'orange',
#        align = 'center',width = np.diff(bins)[0],capsize = 5)        
reslow[reslow<0.1] = 0
ecolow[ecolow<0.1] = 0

plt.errorbar(mhbins+(np.diff(mhbins)[0]/2),eco_sub_agnfrac,fmt = '.-',yerr = [ecolow,ecoup],\
             label = 'ECO SEL Dwarfs (excluding RESOLVE-A)', capsize = 5)
plt.errorbar(mhbins+(np.diff(mhbins)[0]/2),res_sub_agnfrac,fmt = '.-',yerr = [reslow,resup],\
             label = 'RESOLVE-A SEL Dwarfs', capsize = 5)

plt.legend(loc = 'upper left')
#plt.xlabel('Normalized environmental density smoothed at 1Mpc scale', fontsize = 15)
plt.xlabel('log(M_bary)', fontsize = 15)
plt.ylabel('Percentage of Dwarf SEL AGN', fontsize = 15)
plt.xticks(mhbins)
#plt.ylim(0.,60)
#plt.xlim(-1, 5)
#
##Density vs. halo Mass
#plt.figure()
##plt.scatter(np.log10(10**eco.logmstar[ecodwarfndx]+10**eco.logmgas[ecodwarfndx]),ecoden)
##plt.scatter(resamh,resaden)
##plt.yscale('log')
#
#bins = np.arange(10.5,15.5,0.5)
#ecoden_med = []
#resaden_med = []
#ecoden_low = []
#ecoden_high = []
#resaden_low = []
#resaden_high = []
#ecoden_std = []
#resaden_std = []
#for i in range(len(bins)-1):
#    ecoden_sub = ecoden[((ecomh < bins[i+1]) & (ecomh > bins[i]))]
#    resaden_sub = resaden[((resamh < bins[i+1]) & (resamh > bins[i]))]
#
#    ecoden_med.append(np.median(ecoden_sub))
#    resaden_med.append(np.median(resaden_sub))
#
#    ecoden_low.append(np.min(ecoden_sub))
#    ecoden_high.append(np.max(ecoden_sub))
#
#    resaden_low.append(np.min(resaden_sub))
#    resaden_high.append(np.max(resaden_sub))
#
#    ecoden_std.append(np.std(ecoden_sub))    
#    resaden_std.append(np.std(resaden_sub))    
##plt.errorbar(bins[:-1],ecoden_med, fmt = 'o', yerr = [ecoden_low,ecoden_high],
##             label = 'ECO SEL Dwarfs (excluding RESOLVE-A)',
##             capsize = 5)
##plt.errorbar(bins[:-1],resaden_med, fmt = 'o', yerr = [resaden_low, resaden_high],
##             label = 'RESOLVE-A SEL Dwarfs', capsize = 5)
#plt.errorbar(bins[:-1],ecoden_med, fmt = '.-', yerr = ecoden_std,
#             label = 'ECO SEL Dwarfs (excluding RESOLVE-A)',
#             capsize = 5)
#plt.errorbar(bins[:-1],resaden_med, fmt = '.-', yerr = resaden_std,
#             label = 'RESOLVE-A SEL Dwarfs', capsize = 5)
#plt.legend()


ecoagn = ecoflag.sftoagn |  ecoflag.defagn | \
                    ecoflag.composite | ecoflag.agntosf
ecoden = eco.den1mpc[ecoagn & (eco.den1mpc > -99) & (eco.logmstar < 9.5)]
ecomh = eco.logmh[ecoagn & (eco.den1mpc > -99) & (eco.logmstar < 9.5)]


resaagn = ecofullflag[~econdx].sftoagn |  ecofullflag[~econdx].defagn | \
                    ecofullflag[~econdx].composite | ecofullflag[~econdx].agntosf
resaden = ecofull.den1mpc[~econdx][resaagn & (ecofull[~econdx].logmstar < 9.5) \
                         & (ecofull[~econdx].den1mpc > -99)]
resamh = ecofull.logmh[~econdx][resaagn & (ecofull[~econdx].logmstar < 9.5) \
                         & (ecofull[~econdx].den1mpc > -99)]
bins = np.arange(10.5,15.5,0.5)
#ecoden_med = []
#resaden_med = []
#ecoden_low = []
#ecoden_high = []
#resaden_low = []
#resaden_high = []
#ecoden_std = []
#resaden_std = []
#for i in range(len(bins)-1):
#    ecoden_sub = ecoden[((ecomh < bins[i+1]) & (ecomh > bins[i]))]
#    resaden_sub = resaden[((resamh < bins[i+1]) & (resamh > bins[i]))]
#
#    ecoden_med.append(np.median(ecoden_sub))
#    resaden_med.append(np.median(resaden_sub))
#
#    ecoden_low.append(np.min(ecoden_sub))
#    ecoden_high.append(np.max(ecoden_sub))
#
#    resaden_low.append(np.min(resaden_sub))
#    resaden_high.append(np.max(resaden_sub))
#
#    ecoden_std.append(np.std(ecoden_sub))    
#    resaden_std.append(np.std(resaden_sub))    
##plt.errorbar(bins[:-1],ecoden_med, fmt = 'o', yerr = [ecoden_low,ecoden_high],
##             label = 'ECO SEL Dwarfs (excluding RESOLVE-A)',
##             capsize = 5)
##plt.errorbar(bins[:-1],resaden_med, fmt = 'o', yerr = [resaden_low, resaden_high],
##             label = 'RESOLVE-A SEL Dwarfs', capsize = 5)
#plt.errorbar(bins[:-1],ecoden_med, fmt = '.-', yerr = ecoden_std,
#             label = 'ECO SEL Dwarf AGN (excluding RESOLVE-A)',
#             capsize = 5)
#plt.errorbar(bins[:-1],resaden_med, fmt = '.-', yerr = resaden_std,
#             label = 'RESOLVE-A SEL Dwarf AGN', capsize = 5)
#plt.legend()



resnew = pd.read_csv('RESOLVE_liveMay2020.csv')
resnew.index = resnew.name
res['dens_s'] = resnew.loc[res.name]['dens_s']
resmh = res.logmh[(res.logmstar < 9.5) & (res['dens_s'] > 0)]
res_lssdens = res['dens_s'][(res.logmstar < 9.5) & (res['dens_s'] > 0)]
#plt.figure()
#plt.hist(res_lssdens, bins = 'fd')
#

#eco['dens_s'] = np.zeros(len(eco))
#eco['dens_s'].loc[list(ecocat.econames[dat_ndx])]= np.log10(eco_lssdens_dat.dm_dens[dat_ndx])
#eco_lssdens = eco['dens_s'][(eco.logmstar < 9.5) & (eco['dens_s'] > 0)]
#ecomh = eco['logmh'][(eco.logmstar < 9.5) & (eco['dens_s'] > 0)]
#
#plt.figure()
##plt.scatter(np.log10(10**eco.logmstar[ecodwarfndx]+10**eco.logmgas[ecodwarfndx]),ecoden)
##plt.scatter(resamh,resaden)
##plt.yscale('log')
#
#bins = np.arange(10.5,15.5,0.5)
#ecoden_med = []
#resaden_med = []
#ecoden_low = []
#ecoden_high = []
#resaden_low = []
#resaden_high = []
#ecoden_std = []
#resaden_std = []
#for i in range(len(bins)-1):
#    ecoden_sub = eco_lssdens[((ecomh < bins[i+1]) & (ecomh > bins[i]))]
#    resaden_sub = res_lssdens[((resmh < bins[i+1]) & (resmh > bins[i]))]
#
#    ecoden_med.append(np.median(ecoden_sub))
#    resaden_med.append(np.median(resaden_sub))
#
#    ecoden_low.append(np.min(ecoden_sub))
#    ecoden_high.append(np.max(ecoden_sub))
#
#    resaden_low.append(np.min(resaden_sub))
#    resaden_high.append(np.max(resaden_sub))
#
#    ecoden_std.append(np.std(ecoden_sub))    
#    resaden_std.append(np.std(resaden_sub))    
##plt.errorbar(bins[:-1],ecoden_med, fmt = 'o', yerr = [ecoden_low,ecoden_high],
##             label = 'ECO SEL Dwarfs (excluding RESOLVE-A)',
##             capsize = 5)
##plt.errorbar(bins[:-1],resaden_med, fmt = 'o', yerr = [resaden_low, resaden_high],
##             label = 'RESOLVE-A SEL Dwarfs', capsize = 5)
#plt.errorbar(bins[:-1],ecoden_med, fmt = '.-', yerr = ecoden_std,
#             label = 'ECO SEL Dwarfs (excluding RESOLVE-A)',
#             capsize = 5)
#plt.errorbar(bins[:-1],resaden_med, fmt = '.-', yerr = resaden_std,
#             label = 'RESOLVE-A SEL Dwarfs', capsize = 5)
#plt.xlabel('Log(Halo Mass)')
#plt.ylabel('Typical Large Scale Structure Density')
#plt.legend()


#ecoagn = ecoflag.sftoagn |  ecoflag.defagn | \
#                    ecoflag.composite | ecoflag.agntosf
#eco_lssdens = eco['dens_s'][ecoagn & (eco.logmstar < 9.5) & (eco['dens_s'] > 0)]
#ecomh = eco['logmh'][ecoagn & (eco.logmstar < 9.5) & (eco['dens_s'] > 0)]
#
#
#resagn = resflag.sftoagn |  resflag.defagn | \
#                    resflag.composite | resflag.agntosf
#resmh = res.logmh[resagn & (res.logmstar < 9.5) \
#                         & (res.dens_s > 0)]
#res_lssdens = res.dens_s[resagn & (res.logmstar < 9.5) \
#                         & (res.dens_s > 0)]
#bins = np.arange(10.5,15.5,0.5)
#ecoden_med = []
#resaden_med = []
#ecoden_low = []
#ecoden_high = []
#resaden_low = []
#resaden_high = []
#ecoden_std = []
#resaden_std = []
#for i in range(len(bins)-1):
#    ecoden_sub = eco_lssdens[((ecomh < bins[i+1]) & (ecomh > bins[i]))]
#    resaden_sub = res_lssdens[((resmh < bins[i+1]) & (resmh > bins[i]))]
#
#    ecoden_med.append(np.median(ecoden_sub))
#    resaden_med.append(np.median(resaden_sub))
#
#    ecoden_low.append(np.min(ecoden_sub))
#    ecoden_high.append(np.max(ecoden_sub))
#
#    resaden_low.append(np.min(resaden_sub))
#    resaden_high.append(np.max(resaden_sub))
#
#    ecoden_std.append(np.std(ecoden_sub))    
#    resaden_std.append(np.std(resaden_sub))    
##plt.errorbar(bins[:-1],ecoden_med, fmt = 'o', yerr = [ecoden_low,ecoden_high],
##             label = 'ECO SEL Dwarfs (excluding RESOLVE-A)',
##             capsize = 5)
##plt.errorbar(bins[:-1],resaden_med, fmt = 'o', yerr = [resaden_low, resaden_high],
##             label = 'RESOLVE-A SEL Dwarfs', capsize = 5)
#plt.errorbar(bins[:-1],ecoden_med, fmt = '.-', yerr = ecoden_std,
#             label = 'ECO SEL Dwarf AGN (excluding RESOLVE-A)',
#             capsize = 5)
#plt.errorbar(bins[:-1],resaden_med, fmt = '.-', yerr = resaden_std,
#             label = 'RESOLVE-A SEL Dwarf AGN', capsize = 5)
#plt.legend()

#plt.figure()
##c = ['g','b',r','m','k']
#res_sub_agnfrac, reslow, resup, \
#        eco_sub_agnfrac, ecolow, ecoup,mhbins = \
#        agnfrac_dens('dwarf', 'dens_s')
#reslow[reslow<0.1] = 0
#ecolow[ecolow<0.1] = 0
#reslow[reslow<0.1] = 0
#ecolow[ecolow<0.1] = 0
#
#plt.errorbar(mhbins+0.25,eco_sub_agnfrac,fmt = '.-',yerr = [ecolow,ecoup],\
#             label = 'ECO SEL Dwarfs (excluding RESOLVE-A)', capsize = 5)
#plt.errorbar(mhbins+0.25,res_sub_agnfrac,fmt = '.-',yerr = [reslow,resup],\
#             label = 'RESOLVE-A SEL Dwarfs', capsize = 5)
#
#plt.legend(loc = 'upper left')
#plt.xlabel('Large Scale Structure Density', fontsize = 15)
#plt.ylabel('Percentage of Dwarf SEL AGN', fontsize = 15)
##plt.ylim(0.,60)
#plt.xlim(9, 14)

eco_mhgrp = []
for i in np.unique(eco.grp):
    eco_mhgrp.append((i, eco.logmh[(eco.grp == i)][0]))
ecoparent_mhgrp = []
for i in np.unique(ecoparent.grp):
    ecoparent_mhgrp.append((i, ecoparent.logmh[(ecoparent.grp == i)][0]))
res_mhgrp = []
for i in np.unique(res.grp):
    res_mhgrp.append((i, res.logmh[(res.grp == i)][0]))
resparent_mhgrp = []
for i in np.unique(resparent.grp):
    resparent_mhgrp.append((i, resparent.logmh[(resparent.grp == i)][0]))

ecoagn_mhgrp = []
for i in np.unique(eco.grp):
    ecoagn_mhgrp.append((i, eco.logmh[(eco.grp == i)][0]))
resagn_mhgrp = []
for i in np.unique(res.grp):
    res_mhgrp.append((i, res.logmh[(res.grp == i)][0]))

#plt.figure()
#for i in np.unique(res.grp):
#    grpno = i #res.grp[i]
#    agningroup = (resdwarfagn.grp == grpno)
#    selingroup = (res.grp == grpno)
#    percent = 100.0*(np.sum(agningroup))/np.sum(selingroup)
#    low, up = 100*binom_conf_interval(np.sum(agningroup),\
#                                                  np.sum(selingroup))
#    low_err = percent - low
#    up_err = up - percent
#    ndx = np.where(res.grp == i)[0][0]
#    plt.errorbar(res.logmh[ndx],percent,fmt = 'bo', alpha = 0.3,
#                 yerr = [[low_err], [up_err]])
#for i in np.unique(eco.grp):# range(len(eco)): [454]:#agngrps: #
#    grpno = i #eco.grp[i]
#    agningroup = (ecodwarfagn.grp == grpno)
#    selingroup = (eco.grp == grpno)
#    percent = 100.0*(np.sum(agningroup))/np.sum(selingroup)
#    low, up = 100*binom_conf_interval(np.sum(agningroup),\
#                                                  np.sum(selingroup))
#    low_err = percent - low
#    up_err = up - percent
#    ndx = np.where(eco.grp == i)[0][0]
#    plt.errorbar(eco.logmh[ndx],percent,fmt = 'ro', alpha = 0.3,
#                 yerr = [[low_err], [up_err]])
#nonsels = [x for x in np.unique(ecoparent.grp) if x not in np.unique(eco.grp)]
#for i in nonsels:# range(len(eco)):
#    grpno = i #eco.grp[i]
#    percent = 0.0
##    low, up = 100*binom_conf_interval(np.sum(agningroup),\
##                                                  np.sum(selingroup))
##    low_err = percent - low
##    up_err = up - percent
#    ndx = np.where(ecoparent.grp == i)[0][0]
#    plt.plot(ecoparent.logmh[ndx],percent,'kx')
#plt.xlabel('Halo Mass')
#plt.ylabel('% of SEL AGN in each group')


agngrps = np.array(ecodwarfagn.grp[(ecodwarfagn.logmh > 13.0) & \
                                    (ecodwarfagn.logmh < 13.5)])
agngrps_names = np.array(ecodwarfagn.name[(ecodwarfagn.logmh > 13.0) & \
                                    (ecodwarfagn.logmh < 13.5)])
fig = plt.figure()
ax = fig.add_subplot(111)
#ax.plot(ras, decs, 'ko', ms=3, mfc = 'none', label='Galaxy')
ax.set_xlabel('RA (hours)')
ax.set_ylabel('Dec (degrees)')
ax.set_aspect(aspect=0.0667)
ax.plot(ecoparent.radeg/15,ecoparent.dedeg,'k.', alpha  = 0.1,mew = 0)
ax.plot(ecofull.radeg/15,ecofull.dedeg,'k.', alpha  = 0.3,mew = 0)
agngrps = ecodwarfagn.grp[ecodwarfagn.logmbary < 9.2]
#for i in agngrps:
#    gals = (ecoparent.grp == i)
#    #plt.figure('ECO Group '+str(i))
#    plt.plot(ecoparent.radeg[gals]/15,ecoparent.dedeg[gals],'bo')
#    central = (ecoparent[gals].fc == 1)
#    #plt.plot(ecoparent[gals].radeg[central]/15,ecoparent[gals].dedeg[central], \
#    #         'bo',ms = 10)
#    
#    dwarfagningrp = (ecodwarfagn.grp == i)
#    #plt.plot(ecodwarfagn.radeg[dwarfagningrp]/15,ecodwarfagn.dedeg[dwarfagningrp],'ro')
    
plt.xlabel('RA (hours)')
plt.ylabel('Dec')
plt.plot(ecosatdwarf.radeg/15,ecosatdwarf.dedeg,'bo')
plt.plot(ecosatdwarfagn.radeg/15,ecosatdwarfagn.dedeg,'ro')
plt.plot(ressatdwarf.radeg/15,ressatdwarf.dedeg,'go')
plt.plot(ressatdwarfagn.radeg/15,ressatdwarfagn.dedeg,'mo')
xaxis = np.arange(8.5,16.5,0.1)
ax.plot(xaxis, 5*np.ones(len(xaxis)),'r-.')

normgrps = np.array(eco.grp[(eco.logmh > 13.0) & \
                                    (eco.logmh < 13.5)])
normgrps = [x for x in normgrps if x not in agngrps]
#for i in agngrps: #normgrps[:4]:
#    gals = (ecoparent.grp == i)
#    plt.figure()
#    plt.title('ECO Group '+str(i))
#    plt.plot(ecoparent.radeg[gals],ecoparent.dedeg[gals],'ko')
#    central = (ecoparent[gals].fc == 1)
#    plt.plot(ecoparent[gals].radeg[central],ecoparent[gals].dedeg[central], \
#             'ko', mfc = 'none', ms = 10)
#    
#    sels = (eco.grp == i)
#    plt.plot(eco.radeg[sels],eco.dedeg[sels],'bo')
#    
#    dwarfagningrp = (ecodwarfagn.grp == i)
#    plt.plot(ecodwarfagn.radeg[dwarfagningrp],ecodwarfagn.dedeg[dwarfagningrp],'ro')
#    
#    agningrp = (ecogiantagn.grp == i)
#    plt.plot(ecogiantagn.radeg[agningrp],ecogiantagn.dedeg[agningrp],'go')
#    
#    
#    plt.xlabel('RA')
#    plt.ylabel('Dec')
#
resparent_mhgrp = np.array(resparent_mhgrp)
ecoparent_mhgrp = np.array(ecoparent_mhgrp)
res_mhgrp = np.array(res_mhgrp)
eco_mhgrp = np.array(eco_mhgrp)

bins = np.arange(10.5,15.5,0.5)
fig = plt.figure()
ax1 = plt.subplot(121)
ax1.hist(ecoparent_mhgrp[:,1], bins = bins, histtype = 'step', 
         edgecolor = 'k', lw = 3, label = 'ECO Parent Sample (excluding RESOLVE-A)')#,normed = True)
ax1.hist(eco_mhgrp[:,1], bins = bins, histtype = 'step', hatch = 'x',
         edgecolor = 'orange', lw = 2, label = 'ECO SELs (excluding RESOLVE-A)')#, normed = True)
ax1.legend()
ax1.set_xlabel('log(M$_{halo}$/M$_{\odot}$)', fontsize = 15) 
ax1.set_ylabel('Number of groups', fontsize = 15)
ax2 = plt.subplot(122)
ax2.hist(resparent_mhgrp[:,1], bins = bins, histtype = 'step',  
         edgecolor = 'black', lw = 3, label = 'RESOLVE Parent Sample')#, normed = True)
ax2.hist(res_mhgrp[:,1], bins = bins, histtype = 'step', hatch = 'x',
         edgecolor = 'orange', lw = 2, label = 'RESOLVE SELs')#, normed = True)
ax2.legend()
ax2.set_xlabel('log(M$_{halo}$/M$_{\odot}$)', fontsize = 15)
ax2.set_ylabel('Number of groups', fontsize = 15)

ecog3full = pd.read_csv('ecog3.csv')
ecog3full.index = ecog3full.name
ecog3ndx = (((ecog3full.cz < 4500) | (ecog3full.dedeg < 0) | (ecog3full.dedeg > 5)) | \
         (ecog3full.resname == 'notinresolve'))
ecog3 = ecog3full[ecog3ndx]

agninlargegrp = list(ecodwarfagn.name[(ecodwarfagn.logmh > 13.0) & \
                                    (ecodwarfagn.logmh < 13.5)])

###############################################################################
bins = np.arange(10.5,15,0.5)
plt.figure()
plt.hist(ecosatdwarf.logmh, bins = bins, histtype = 'step', lw = 3, 
         color = 'blue',  normed = True, label = 'ECO Dwarf Sat. (excluding RES-A)')
plt.hist(ecosatdwarfagn.logmh, bins = bins, histtype = 'step', lw = 3, 
         normed = True, color = 'red', label = 'ECO Dwarf Sat. AGN')
plt.hist(ressatdwarf.logmh, bins = bins, histtype = 'step', lw = 3, 
         color = 'green',  normed = True, label = 'RESOLVE Dwarf Sat.')
plt.legend()

resadwarfdenndx = ~econdx & (ecofull.logmstar < 9.5) & (ecofull.den1mpc > -99)
resaden = ecofull.den1mpc[resadwarfdenndx]

###############################################################################
bins = np.arange(0,4,0.5)
plt.figure()
plt.hist(ecosatdwarf.den1mpc, bins = bins, histtype = 'step', lw = 3, 
         color = 'blue',  normed = True, label = 'ECO Dwarf Sat. (excluding RES-A)')
plt.hist(ecosatdwarfagn.den1mpc, bins = bins, histtype = 'step', lw = 3, 
         normed = True, color = 'red', label = 'ECO Dwarf Sat. AGN')
plt.hist(resaden, bins = bins, histtype = 'step', lw = 3, 
         color = 'green',  normed = True , label = 'RESOLVE Dwarf Sat.')
plt.legend()

###############################################################################

ecohasnr = ecosatdwarf.h_alpha_flux/ecosatdwarf.h_alpha_flux_err
ecooisnr = ecosatdwarf.oi_6300_flux/ecosatdwarf.oi_6300_flux_err

ecoagnhasnr = ecosatdwarfagn.h_alpha_flux/ecosatdwarfagn.h_alpha_flux_err
ecoagnoisnr = ecosatdwarfagn.oi_6300_flux/ecosatdwarfagn.oi_6300_flux_err

reshasnr = ressatdwarf.h_alpha_flux/ressatdwarf.h_alpha_flux_err
resoisnr = ressatdwarf.oi_6300_flux/ressatdwarf.oi_6300_flux_err

plt.figure()
plt.hist(ecohasnr, bins = bins, histtype = 'step', lw = 3, 
         color = 'blue',  normed = True, label = 'ECO Dwarf Sat. (excluding RES-A)')
plt.hist(ecoagnhasnr, bins = bins, histtype = 'step', lw = 3, 
         normed = True, color = 'red', label = 'ECO Dwarf Sat. AGN')
plt.hist(reshasnr, bins = bins, histtype = 'step', lw = 3, 
         color = 'green',  normed = True , label = 'RESOLVE Dwarf Sat.')
plt.legend()

plt.figure()
plt.hist(ecooisnr, bins = bins, histtype = 'step', lw = 3, 
         color = 'blue',  normed = True, label = 'ECO Dwarf Sat. (excluding RES-A)')
plt.hist(ecoagnoisnr, bins = bins, histtype = 'step', lw = 3, 
         normed = True, color = 'red', label = 'ECO Dwarf Sat. AGN')
plt.hist(resoisnr, bins = bins, histtype = 'step', lw = 3, 
         color = 'green',  normed = True , label = 'RESOLVE Dwarf Sat.')
plt.legend()
    
###############################################################################
#Density vs. Baryonic Mass
###############################################################################
#ecodenndx = (eco.logmstar < 9.5) & (eco.den1mpc > -99)
#ecoden = eco.den1mpc[ecodenndx]
#ecobary = np.log10(10**eco.logmstar[ecodenndx] + 10**eco.logmgas[ecodenndx])
#resadenndx = ~econdx & (ecofull.logmstar < 9.5) & (ecofull.den1mpc > -99)
#resaden = ecofull.den1mpc[resadenndx]
#resabary = np.log10(10**ecofull.logmstar + 10**ecofull.logmgas)[resadenndx]
#
#
#plt.figure()
##plt.scatter(np.log10(10**eco.logmstar[ecodwarfndx]+10**eco.logmgas[ecodwarfndx]),ecoden)
##plt.scatter(resamh,resaden)
##plt.yscale('log')
#
#bins = np.arange(8.5,12.0,0.5)
#ecoden_med = []
#resaden_med = []
#ecoden_low = []
#ecoden_high = []
#resaden_low = []
#resaden_high = []
#ecoden_std = []
#resaden_std = []
#for i in range(len(bins)-1):
#    ecoden_sub = ecoden[((ecobary< bins[i+1]) & (ecobary> bins[i]))]
#    resaden_sub = resaden[((resabary< bins[i+1]) & (resabary> bins[i]))]
#
#    ecoden_med.append(np.median(ecoden_sub))
#    resaden_med.append(np.median(resaden_sub))
#
#    ecoden_low.append(np.min(ecoden_sub))
#    ecoden_high.append(np.max(ecoden_sub))
#
#    resaden_low.append(np.min(resaden_sub))
#    resaden_high.append(np.max(resaden_sub))
#
#    ecoden_std.append(np.std(ecoden_sub))    
#    resaden_std.append(np.std(resaden_sub))    
##plt.errorbar(bins[:-1],ecoden_med, fmt = 'o', yerr = [ecoden_low,ecoden_high],
##             label = 'ECO SEL Dwarfs (excluding RESOLVE-A)',
##             capsize = 5)
##plt.errorbar(bins[:-1],resaden_med, fmt = 'o', yerr = [resaden_low, resaden_high],
##             label = 'RESOLVE-A SEL Dwarfs', capsize = 5)
#plt.errorbar(bins[:-1],ecoden_med, fmt = '.-', yerr = ecoden_std,
#             label = 'ECO SEL Dwarfs (excluding RESOLVE-A)',
#             capsize = 5)
#plt.errorbar(bins[:-1],resaden_med, fmt = '.-', yerr = resaden_std,
#             label = 'RESOLVE-A SEL Dwarfs', capsize = 5)
#plt.legend()
#
#
#ecoagn = ecoflag.sftoagn |  ecoflag.defagn | \
#                    ecoflag.composite | ecoflag.agntosf
#
#ecodenndx = (ecoagn) & (eco.logmstar < 9.5) & (eco.den1mpc > -99)
#ecoden = eco.den1mpc[ecodenndx]
#ecobary = np.log10(10**eco.logmstar[ecodenndx] + 10**eco.logmgas[ecodenndx])
#resaagn = ecofullflag[~econdx].sftoagn |  ecofullflag[~econdx].defagn | \
#                    ecofullflag[~econdx].composite | ecofullflag[~econdx].agntosf
#resadenndx = resaagn & (ecofull[~econdx].logmstar < 9.5) \
#                         & (ecofull[~econdx].den1mpc > -99)
#resaden = ecofull.den1mpc[~econdx][resadenndx]
#resabary = np.log10(10**ecofull.logmstar[~econdx] + \
#                    10**ecofull.logmgas[~econdx])[resadenndx]
#
#
#ecoden_med = []
#resaden_med = []
#ecoden_low = []
#ecoden_high = []
#resaden_low = []
#resaden_high = []
#ecoden_std = []
#resaden_std = []
#for i in range(len(bins)-1):
#    ecoden_sub = ecoden[((ecobary< bins[i+1]) & (ecobary> bins[i]))]
#    resaden_sub = resaden[((resabary< bins[i+1]) & (resabary> bins[i]))]
#
#    ecoden_med.append(np.median(ecoden_sub))
#    resaden_med.append(np.median(resaden_sub))
#
#    ecoden_low.append(np.min(ecoden_sub))
#    ecoden_high.append(np.max(ecoden_sub))
#
#    resaden_low.append(np.min(resaden_sub))
#    resaden_high.append(np.max(resaden_sub))
#
#    ecoden_std.append(np.std(ecoden_sub))    
#    resaden_std.append(np.std(resaden_sub))    
##plt.errorbar(bins[:-1],ecoden_med, fmt = 'o', yerr = [ecoden_low,ecoden_high],
##             label = 'ECO SEL Dwarfs (excluding RESOLVE-A)',
##             capsize = 5)
##plt.errorbar(bins[:-1],resaden_med, fmt = 'o', yerr = [resaden_low, resaden_high],
##             label = 'RESOLVE-A SEL Dwarfs', capsize = 5)
#plt.errorbar(bins[:-1],ecoden_med, fmt = '.-', yerr = ecoden_std,
#             label = 'ECO SEL Dwarf AGN (excluding RESOLVE-A)',
#             capsize = 5)
#plt.errorbar(bins[:-1],resaden_med, fmt = '.-', yerr = resaden_std,
#             label = 'RESOLVE-A SEL Dwarf AGN', capsize = 5)
#plt.legend()
    
###############################################################################    
# AGN % as a function of Halo mass, binned by Declination
###############################################################################    
#def subset_agnfrac(lowdec,highdec,galtype = 'all'):
#    if galtype == 'all':
#        res_sub_ndx = (res.dedeg > lowdec) & (res.dedeg < highdec)
#        eco_sub_ndx = (eco.dedeg > lowdec) & (eco.dedeg < highdec)
#    elif galtype == 'dwarf':
#        res_sub_ndx = (res.dedeg > lowdec) & (res.dedeg < highdec) & (res.logmstar < 9.5)
#        eco_sub_ndx = (eco.dedeg > lowdec) & (eco.dedeg < highdec) & (eco.logmstar < 9.5)
#    
#    res_sub = res[res_sub_ndx]
#    res_sub_flag = resflag[res_sub_ndx]
#    eco_sub = eco[eco_sub_ndx]
#    eco_sub_flag = ecoflag[eco_sub_ndx]
#    
#    res_sub_agn = res_sub_flag.sftoagn |  res_sub_flag.defagn | \
#                    res_sub_flag.composite | res_sub_flag.agntosf
#    eco_sub_agn = eco_sub_flag.sftoagn |  eco_sub_flag.defagn | \
#                    eco_sub_flag.composite | eco_sub_flag.agntosf
#
#    res_sub_agnfrac = 100.0*np.sum(res_sub_agn)/len(res_sub)
#    eco_sub_agnfrac = 100.0*np.sum(eco_sub_agn)/len(eco_sub)
#
#    if len(res_sub) > 0:
#        reslowlim, resuplim = 100*binom_conf_interval(np.sum(res_sub_agn),len(res_sub))
#        resup= resuplim - res_sub_agnfrac
#        reslow = res_sub_agnfrac - reslowlim 
#    else:
#        resup= 0
#        reslow = 0
#    ecolowlim, ecouplim = 100*binom_conf_interval(np.sum(eco_sub_agn),len(eco_sub))
#    ecoup= ecouplim - eco_sub_agnfrac
#    ecolow = eco_sub_agnfrac - ecolowlim 
#
#    return res_sub_agnfrac, reslow, resup, \
#            eco_sub_agnfrac, ecolow, ecoup
#
##Dec -1 to 5
#res_sub_agnfrac = np.zeros(9)
#lowerr_res_sub_agnfrac = np.zeros((9))
#uperr_res_sub_agnfrac = np.zeros((9))
#eco_sub_agnfrac = np.zeros(9)
#lowerr_eco_sub_agnfrac = np.zeros((9))
#uperr_eco_sub_agnfrac = np.zeros((9))
#lowdec = np.zeros(9)
#highdec = np.zeros(9)
#for i in range(9):
#    if i == 0:
#        lowdec[i] = -1
#    else:
#        lowdec[i] = i*5
#    highdec[i] = lowdec[i]+5
#    
#    res_sub_agnfrac[i], lowerr_res_sub_agnfrac[i], uperr_res_sub_agnfrac[i], \
#    eco_sub_agnfrac[i], lowerr_eco_sub_agnfrac[i], uperr_eco_sub_agnfrac[i] = \
#                        subset_agnfrac(lowdec[i], highdec[i])
#
#plt.figure()
#plt.errorbar(highdec,res_sub_agnfrac, fmt = 'bo', yerr = [lowerr_res_sub_agnfrac,\
#                                                   uperr_res_sub_agnfrac])
#plt.errorbar(highdec,eco_sub_agnfrac, fmt = 'ro', yerr = [lowerr_eco_sub_agnfrac,\
#                                                    uperr_eco_sub_agnfrac])
#plt.xlabel('Dec')
#plt.ylabel('AGN Fraction')
#
#res_sub_agnfrac = np.zeros(9)
#err_res_sub_agnfrac = np.zeros((2,9))
#eco_sub_agnfrac = np.zeros(9)
#err_eco_sub_agnfrac = np.zeros((2,9))
#lowdec = np.zeros(9)
#highdec = np.zeros(9)
#for i in range(9):
#    if i == 0:
#        lowdec[i] = -1
#    else:
#        lowdec[i] = i*5
#    highdec[i] = lowdec[i]+5
#    
#    res_sub_agnfrac[i], lowerr_res_sub_agnfrac[i], uperr_res_sub_agnfrac[i], \
#    eco_sub_agnfrac[i], lowerr_eco_sub_agnfrac[i], uperr_eco_sub_agnfrac[i] = \
#                        subset_agnfrac(lowdec[i], highdec[i], 'dwarf')
#
#plt.figure()
#plt.errorbar(highdec,res_sub_agnfrac, fmt = 'bo', yerr = [lowerr_res_sub_agnfrac,\
#                                                   uperr_res_sub_agnfrac])
#plt.errorbar(highdec,eco_sub_agnfrac, fmt = 'ro', yerr = [lowerr_eco_sub_agnfrac,\
#                                                    uperr_eco_sub_agnfrac])
#plt.xlabel('Dec')
#plt.ylabel('Dwarf AGN Fraction')
#def subset_agnfrac_env(lowdec,highdec,galtype = 'all'):
#    if galtype == 'all':
#        res_sub_ndx = (res.dedeg > lowdec) & (res.dedeg < highdec)
#        eco_sub_ndx = (eco.dedeg > lowdec) & (eco.dedeg < highdec)
#    elif galtype == 'dwarf':
#        res_sub_ndx = (res.dedeg > lowdec) & (res.dedeg < highdec) & \
#        (res.logmstar < 9.5)
#        eco_sub_ndx = (eco.dedeg > lowdec) & (eco.dedeg < highdec) & \
#        (eco.logmstar < 9.5)
#    
#    res_sub = res[res_sub_ndx]
#    res_sub_flag = resflag[res_sub_ndx]
#    eco_sub = eco[eco_sub_ndx]
#    eco_sub_flag = ecoflag[eco_sub_ndx]
#    
#    res_sub_agn = res_sub_flag.sftoagn |  res_sub_flag.defagn | \
#                    res_sub_flag.composite | res_sub_flag.agntosf
#    eco_sub_agn = eco_sub_flag.sftoagn |  eco_sub_flag.defagn | \
#                    eco_sub_flag.composite | eco_sub_flag.agntosf
#
#    #res_sub_agnfrac = 100.0*np.sum(res_sub_agn)/len(res_sub)
#    #eco_sub_agnfrac = 100.0*np.sum(eco_sub_agn)/len(eco_sub)
#
#    
#    res_sub_env = np.zeros(5)
#    eco_sub_env = np.zeros(5)
#    res_sub_agn_env = np.zeros(5)
#    eco_sub_agn_env = np.zeros(5)
#    mh_bins = np.zeros(5)
#    for i in range(5):
#        lowmh = 10.5 + i*1
#        highmh = lowmh + 1
#        
#        mh_bins[i] = highmh
#        res_env_ndx = (res_sub.logmh > lowmh) & (res_sub.logmh < highmh)
#        res_sub_env[i] = len(res_sub[res_env_ndx])
#        res_sub_agn_env[i] = np.sum(res_sub_agn & res_env_ndx)
#        
#        eco_env_ndx = (eco_sub.logmh > lowmh) & (eco_sub.logmh < highmh)
#        eco_sub_env[i] = len(eco_sub[eco_env_ndx])
#        eco_sub_agn_env[i] = np.sum(eco_sub_agn & eco_env_ndx)
#    
#    
#    res_sub_agnfrac = 100.0*(res_sub_agn_env)/(res_sub_env)
#    eco_sub_agnfrac = 100.0*(eco_sub_agn_env)/(eco_sub_env)
#
#    reslow = np.zeros(5)
#    resup = np.zeros(5)
#    ecolow = np.zeros(5)
#    ecoup = np.zeros(5)
#    for i in range(len(res_sub_agn_env)):
#        if res_sub_agn_env[i] > 0:
#        
#            reslowlim, resuplim = 100*binom_conf_interval(res_sub_agn_env[i],\
#                                                      res_sub_env[i])
#            resup[i]= resuplim - res_sub_agnfrac[i]
#            reslow[i] = res_sub_agnfrac[i] - reslowlim 
#        else:
#            resup[i]= 0
#            reslow[i] = 0
#    for i in range(len(res_sub_agn_env)):
#        if eco_sub_agn_env[i] > 0:
#
#            ecolowlim, ecouplim = 100*binom_conf_interval(eco_sub_agn_env[i],\
#                                                          eco_sub_env[i])
#            ecoup[i]= ecouplim - eco_sub_agnfrac[i]
#            ecolow[i] = eco_sub_agnfrac[i] - ecolowlim 
#        else:
#            ecoup[i]= 0
#            ecolow[i] = 0
#    
#    res_sub_agnfrac = res_sub_agn_env
#    eco_sub_agnfrac = eco_sub_agn_env
#    ecoup= ecoup/100*len(eco_sub)
#    ecolow = ecolow/100*len(eco_sub)
#    resup= ecoup/100*len(res_sub)
#    reslow = ecolow/100*len(res_sub)
#    
#    
#    return res_sub_agnfrac, reslow, resup, \
#            eco_sub_agnfrac, ecolow, ecoup,mh_bins 
#plt.figure()
##c = ['g','b',r','m','k']
#lowdec = [-1,5,10,20,30,40]
#highdec = [5,10,20,30,40,50]
#for i in range(len(lowdec)):
##    if i == 0:
##        lowdec[i] = -1
##    else:
##        lowdec[i] = i*5
##    highdec[i] = lowdec[i]+5
#    
#    res_sub_agnfrac, reslow, resup, \
#            eco_sub_agnfrac, ecolow, ecoup,mhbins = \
#            subset_agnfrac_env(lowdec[i], highdec[i], 'dwarf')
#            
#    plt.errorbar(mhbins-1,eco_sub_agnfrac,fmt = 'o',yerr = [ecolow,ecoup],\
#                 label = 'ECO Dec '+str(lowdec[i])+' to '+str(highdec[i]))
#    if i == 0:
#        plt.errorbar(mhbins-1.1,res_sub_agnfrac,fmt = 'o',yerr = [reslow,resup],\
#                 label = 'RESOLVE Dec '+str(lowdec[i])+' to '+str(highdec[i]))
#plt.legend(loc = 'upper left')
#plt.xlabel('Log(Halo Mass)')
#plt.ylabel('Dwarf AGN Fraction')
##plt.ylim(0.,60)
#plt.xlim(10.0,15.0)