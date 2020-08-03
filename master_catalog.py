
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 18:59:10 2019

@author: mugdhapolimera

Making a master catalog from MPA-JHU, Portsmouth, NSA catalogs
"""

import pandas as pd
import os
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.io.idl import readsav
import numpy as np
from astropy.stats import binom_conf_interval

os.chdir('C:\Users\mugdhapolimera\github\SDSS_spectra/')
resolveflag = 1
ecoflag = 0
fullflag = 0
makemaster = 0
if resolveflag: 
    resolve = pd.read_csv('RESOLVE_full_raw.csv', index_col = 'name')
    jhuflag = pd.read_csv('resolve_emlineclass_dext_snr5_jhu.csv', index_col = 'galname')
    portflag = pd.read_csv('resolve_emlineclass_dext_snr5_port.csv', index_col = 'galname')
    nsaflag = pd.read_csv('resolve_emlineclass_dext_snr5_nsa.csv', index_col = 'galname')
    
    port = pd.read_csv('RESOLVE_full_snr5_dext_port.csv', index_col = 'name')#[portflag.sftoagn]
    jhu = pd.read_csv('RESOLVE_full_snr5_dext_jhu.csv', index_col = 'name')#[jhuflag.sftoagn]
    nsa = pd.read_csv('RESOLVE_full_snr5_dext_nsa.csv', index_col = 'name')#.loc[nsaflag.index.values[nsaflag.sftoagn]]

#if fullflag: 
#    resolve = pd.read_csv('ECO+RESOLVE_full_blend_dext_new.csv', index_col = 'name')
#    jhuflag = pd.read_csv('resolve_emlineclass_full_snr5_jhu.csv', index_col = 'galname')
#    portflag = pd.read_csv('resolve_emlineclass_full_snr5_port.csv', index_col = 'galname')
#    nsaflag = pd.read_csv('resolve_emlineclass_full_snr5_nsa.csv', index_col = 'galname')
#    port = pd.read_csv('ECO+RESOLVE_full_snr5_port.csv', index_col = 'name')#[portflag.sftoagn]
#    jhu = pd.read_csv('ECO+RESOLVE_full_snr5.csv', index_col = 'name')#[jhuflag.sftoagn]
#    nsa = pd.read_csv('NSA_RESOLVE.csv', index_col = 'resname')#.loc[nsaflag.index.values[nsaflag.sftoagn]]

if ecoflag:
#    resolve = pd.read_csv('ECO_live22Oct2018.csv', index_col = 'name')
#    jhuflag = pd.read_csv('eco_emlineclass_full_bary_jhu_new.csv', index_col = 'galname')
#    portflag = pd.read_csv('eco_emlineclass_full_bary_port_new.csv', index_col = 'galname')
#    nsaflag = pd.read_csv('eco_emlineclass_full_snr5_nsa.csv', index_col = 'galname')
#    
#    port = pd.read_csv('ECO_full_bary_port.csv', index_col = 'name')#[portflag.sftoagn]
#    jhu = pd.read_csv('ECO_full_bary_jhu.csv', index_col = 'name')#[jhuflag.sftoagn]
#    nsa = pd.read_csv('NSA_ECO_snr5.csv', index_col = 'name')#.loc[nsaflag.index.values[nsaflag.sftoagn]]
    resolve = pd.read_csv('ECO_full_raw.csv', index_col = 'name')
    os.chdir('C:\Users\mugdhapolimera\github\SDSS_spectra\ECO\SEL')
    jhuflag = pd.read_csv('eco_emlineclass_dext_snr5_jhu.csv', index_col = 'galname')
    portflag = pd.read_csv('eco_emlineclass_dext_snr5_port.csv', index_col = 'galname')
    nsaflag = pd.read_csv('eco_emlineclass_dext_snr5_nsa.csv', index_col = 'galname')
    
    port = pd.read_csv('ECO_full_snr5_dext_port.csv', index_col = 'name')#[portflag.sftoagn]
    jhu = pd.read_csv('ECO_full_dext_snr5_jhu.csv', index_col = 'name')#[jhuflag.sftoagn]
    nsa = pd.read_csv('ECO_full_snr5_dext_nsa.csv', index_col = 'name')#.loc[nsaflag.index.values[nsaflag.sftoagn]]
jhuname = list(jhu.index)
portname = list(port.index)
nsaname = list(nsa.index)


allunq = np.unique(list(jhu.index) + list(port.index)+ list(nsa.index)) 
#sftoagn = df[ambigsel1 & dwarf][['radeg','dedeg']]

#sfagn = pd.read_csv('uniquesfagn.csv')
unique = np.unique(list(jhuflag.index[jhuflag.sftoagn]) + \
                   list(portflag.index[portflag.sftoagn]))
#list(nsaflag.index[nsaflag.sftoagn]) + \
                   
#print('List of Unique SFing-AGN from JHU, Portsmouth and NSA Catalogs')
#print(len(unique))
#print (jhu.loc[[x for x in jhu.index if x in unique]])
#print (port.loc[[x for x in port.index \
#                 if (x in unique) & (x not in jhu.index)]])
#print (nsa.loc[[x for x in nsa.index \
#                if (x in unique) & (x not in jhu.index) & (x not in port.index)]])

df = jhu
df['source'] = 'jhu'
nsakeys = ['nii_6584_flux', 'nii_6584_flux_err', 
           'sii_6731_flux', 'sii_6731_flux_err',
           'oi_6300_flux', 'oi_6300_flux_err', 
           'oiii_5007_flux', 'oiii_5007_flux_err', 
           'h_beta_flux', 'h_beta_flux_err',
           'h_alpha_flux', 'h_alpha_flux_err']
portkeys = ['Flux_NII_6547', 'Flux_NII_6547_Err',
            'Flux_NII_6583', 'Flux_NII_6583_Err', 
           'Flux_SII_6716', 'Flux_SII_6716_Err',
           'Flux_SII_6730', 'Flux_SII_6730_Err', 
           'Flux_OI_6300', 'Flux_OI_6300_Err', 
           'Flux_OIII_5006', 'Flux_OIII_5006_Err',
           'Flux_Hb_4861', 'Flux_Hb_4861_Err', 
           'Flux_Ha_6562', 'Flux_Ha_6562_Err']
jhukeys = ['nii_6548_flux', 'nii_6548_flux_err',
           'nii_6584_flux', 'nii_6584_flux_err', 
           'sii_6717_flux', 'sii_6717_flux_err',
           'sii_6731_flux', 'sii_6731_flux_err',
           'oi_6300_flux', 'oi_6300_flux_err', 
           'oiii_5007_flux', 'oiii_5007_flux_err', 
           'h_beta_flux', 'h_beta_flux_err',
           'h_alpha_flux', 'h_alpha_flux_err']

jhuagn = jhuflag.sftoagn | jhuflag.composite | jhuflag.defagn | jhuflag.agntosf
nsaagn = nsaflag.sftoagn | nsaflag.composite | nsaflag.defagn | nsaflag.agntosf
portagn = portflag.sftoagn | portflag.composite | portflag.defagn | portflag.agntosf

jhudwarf = jhu.logmstar < 9.5
portdwarf = port.logmstar < 9.5
nsadwarf = nsa.logmstar < 9.5

jhudwarfagn = jhudwarf & jhuagn
portdwarfagn = portdwarf & portagn
nsadwarfagn = nsadwarf & nsaagn

print('JHU and Port', 'JHU and NSA', 'NSA and Port')
jhuandport = list(set(list(jhu.index)) & set(list(port.index)))
jhuandnsa = list(set(list(jhu.index)) & set(list(nsa.index)))
nsaandport = list(set(list(nsa.index)) & set(list(port.index)))
print(len(jhuandport), len(jhuandnsa), len(nsaandport))

jhuandportdwarfs = list(set(list(jhu[jhudwarf].index)) & set(list(port[portdwarf].index)))
jhuandnsadwarfs = list(set(list(jhu[jhudwarf].index)) & set(list(nsa[nsadwarf].index)))
nsaandportdwarfs = list(set(list(nsa[nsadwarf].index)) & set(list(port[portdwarf].index)))
print (len(jhuandportdwarfs), len(jhuandnsadwarfs), len(nsaandportdwarfs))

jhuandportdwarfagn = list(set(list(jhu[jhudwarfagn].index)) & set(list(port[portdwarfagn].index)))
jhuandnsadwarfagn = list(set(list(jhu[jhudwarfagn].index)) & set(list(nsa[nsadwarfagn].index)))
nsaandportdwarfagn = list(set(list(nsa[nsadwarfagn].index)) & set(list(port[portdwarfagn].index)))
print (len(jhuandportdwarfagn), len(jhuandnsadwarfagn), len(nsaandportdwarfagn))

print (len(jhuandportdwarfagn), len(jhuandnsadwarfagn), len(nsaandportdwarfagn))
print (100.0*len(jhuandportdwarfagn)/len(jhuandportdwarfs), \
       100.0*binom_conf_interval(len(jhuandportdwarfagn),len(jhuandportdwarfs)) -\
       100.0*len(jhuandportdwarfagn)/len(jhuandportdwarfs))
print (100.0*len(jhuandnsadwarfagn)/len(jhuandnsadwarfs),\
       100.0*binom_conf_interval(len(jhuandnsadwarfagn),len(jhuandnsadwarfs)) -\
       100.0*len(jhuandnsadwarfagn)/len(jhuandnsadwarfs))
print (100.0*len(nsaandportdwarfagn)/len(nsaandportdwarfs),\
       100.0*binom_conf_interval(len(nsaandportdwarfagn),len(nsaandportdwarfs))-\
       100.0*len(nsaandportdwarfagn)/len(nsaandportdwarfs))

print('JHU or Port', 'JHU or NSA', 'NSA or Port')
jhuorport = list(set(list(jhu.index)) | set(list(port.index)))
jhuornsa = list(set(list(jhu.index)) | set(list(nsa.index)))
nsaorport = list(set(list(nsa.index)) | set(list(port.index)))
print(len(jhuorport), len(jhuornsa), len(nsaorport))

jhuorportdwarfs = list(set(list(jhu[jhudwarf].index)) | set(list(port[portdwarf].index)))
jhuornsadwarfs = list(set(list(jhu[jhudwarf].index)) | set(list(nsa[nsadwarf].index)))
nsaorportdwarfs = list(set(list(nsa[nsadwarf].index)) | set(list(port[portdwarf].index)))
print (len(jhuorportdwarfs), len(jhuornsadwarfs), len(nsaorportdwarfs))

jhuorportdwarfagn = list(set(list(jhu[jhudwarfagn].index)) | set(list(port[portdwarfagn].index)))
jhuornsadwarfagn = list(set(list(jhu[jhudwarfagn].index)) | set(list(nsa[nsadwarfagn].index)))
nsaorportdwarfagn = list(set(list(nsa[nsadwarfagn].index)) | set(list(port[portdwarfagn].index)))
print (len(jhuorportdwarfagn), len(jhuornsadwarfagn), len(nsaorportdwarfagn))
print (100.0*len(jhuorportdwarfagn)/len(jhuorportdwarfs), \
       100.0*binom_conf_interval(len(jhuorportdwarfagn),len(jhuorportdwarfs))-\
       100.0*len(jhuorportdwarfagn)/len(jhuorportdwarfs))
print (100.0*len(jhuornsadwarfagn)/len(jhuornsadwarfs),\
       100.0*binom_conf_interval(len(jhuornsadwarfagn),len(jhuornsadwarfs))-\
       100.0*len(jhuornsadwarfagn)/len(jhuornsadwarfs))
print (100.0*len(nsaorportdwarfagn)/len(nsaorportdwarfs),\
       100.0*binom_conf_interval(len(nsaorportdwarfagn),len(nsaorportdwarfs))-\
       100.0*len(nsaorportdwarfagn)/len(nsaorportdwarfs))


print('JHU or Port or NSA')
jhuorportornsa = list(set(list(jhu.index)) | set(list(port.index)) | set(list(nsa.index)))
print (len(jhuorportornsa))

jhuorportornsadwarfs = list(set(list(jhu[jhudwarf].index)) | set(list(port[portdwarf].index)) | set(list(nsa[nsadwarf].index)))
print (len(jhuorportornsadwarfs))

jhuorportornsadwarfagn = list(set(list(jhu[jhudwarfagn].index)) | set(list(port[portdwarfagn].index)) | set(list(nsa[nsadwarfagn].index)))
print (len(jhuorportornsadwarfagn))
print(100.0*len(jhuorportornsadwarfagn)/len(jhuorportornsadwarfs),\
      100.0*binom_conf_interval(len(jhuorportornsadwarfagn),len(jhuorportornsadwarfs))-\
      100.0*len(jhuorportornsadwarfagn)/len(jhuorportornsadwarfs))

print('JHU and Port and NSA')
jhuandportandnsa = list(set(list(jhu.index)) & set(list(port.index)) & set(list(nsa.index)))
print (len(jhuandportandnsa))

jhuandportandnsadwarfs = list(set(list(jhu[jhudwarf].index)) & set(list(port[portdwarf].index)) & set(list(nsa[nsadwarf].index)))
print (len(jhuandportandnsadwarfs))

jhuandportandnsadwarfagn = list(set(list(jhu[jhudwarfagn].index)) & set(list(port[portdwarfagn].index)) & set(list(nsa[nsadwarfagn].index)))
print (len(jhuandportandnsadwarfagn))
print(100.0*len(jhuandportandnsadwarfagn)/len(jhuandportandnsadwarfs),\
      100.0*binom_conf_interval(len(jhuandportandnsadwarfagn),len(jhuandportandnsadwarfs))-\
      100.0*len(jhuandportandnsadwarfagn)/len(jhuandportandnsadwarfs))
#nsaagn = nsaflag.index.values[nsaagnflag]
#portagn = portflag.index.values[portagnflag]
#nsasf = nsaflag.index.values[nsaflag.defstarform]
#portsf = portflag.index.values[portflag.defstarform]
#jhusf = jhuflag.index.values[jhuflag.defstarform]

if makemaster: 

    for gal in allunq:
        #flag = True
        #if gal in jhuagn.index.values:
        #    flag = (~jhuagn.loc[gal]) # | (gal in portagn) | (gal in nsaagn)
        #if flag:
        if gal not in jhu.index.values:        
            if ((gal in port.index.values)) :
                #& (gal not in nsa.index.values)) | \
            #((gal in portagn) & (gal in nsaagn)): 
                if (gal not in df.index.values):
                    print 'port', gal
                    df = df.append(resolve.loc[gal])
                df['source'][gal] = 'port'
                
                for key in range(len(jhukeys)):
                    df.set_value(gal, jhukeys[key], port.loc[gal][portkeys[key]])
            #elif gal in nsa.index.values:
            #    if (gal not in df.index.values):
            #        print 'nsa', gal
            #        df = df.append(resolve.loc[gal])
            #    df['source'][gal] = 'nsa'
            #    df.loc[gal, nsakeys] = nsa.loc[gal][nsakeys]
            #    df.loc[gal,'sii_6717_flux'] = 0.000001
            #    df.loc[gal,'sii_6717_flux_err'] = 0.000001
            #    df.loc[gal,'nii_6548_flux'] = 0.000001
            #    df.loc[gal,'nii_6548_flux_err'] = 0.000001
    if ecoflag:
        df.to_csv('ECO_snr5_master_hasnr5.csv')
    if resolveflag:
        df.to_csv('RESOLVE_snr5_master_new.csv')
    if fullflag:
        df.to_csv('ECO+RESOLVE_snr5_master_new.csv')

    master = df
'''
if resolveflag : 
    master = pd.read_csv('RESOLVE_snr5_master_bary.csv')
if ecoflag : 
    master = pd.read_csv('ECO_snr5_master_bary.csv')
master.index = master.name
num = len(master)
conf = pd.DataFrame({'name': master.index.values, 
                           'jhu': np.zeros(num), 'port': np.zeros(num), 
                           'confidence_level': np.zeros(num)}) #'nsa': np.zeros(num),
conf.index = conf.name
#conf['nsa'].loc[nsaagn] = 1
#conf['nsa'].loc[nsasf] = -1

conf['port'].loc[portagn] = 1
conf['port'].loc[portsf] = -1

conf['jhu'].loc[jhuagn.index.values[jhuagn]] = 1
conf['jhu'].loc[jhusf] = -1

conf['confidence_level'] = conf['port'] + conf['jhu'] #+conf['nsa'] 

final_sample = master.index.values

finalconf = conf.loc[final_sample]

#print(len(np.where(finalconf.confidence_level == 1.0)[0]), 
#      len(np.where(finalconf.confidence_level == 2.0)[0]), 
#      len(np.where(finalconf.confidence_level == 3.0)[0]))

dwarf = master.logmstar < 9.5
agn = np.unique(list(jhuagn.index.values[jhuagn]) + list(portagn))#+ list(nsaagn) + 
                
dwarfagn = master.logmstar.loc[agn] < 9.5
finaldwarf = finalconf.loc[dwarfagn.index.values[dwarfagn]]
print("Confidence Levels for Dwarf AGN found in Master Catalog from "+
      "NSA, JHU, Portsmouth Catalogs")
conf_1 = np.where(finaldwarf.confidence_level == -1.0)[0]
conf0 = np.where(finaldwarf.confidence_level == 0.0)[0]
conf1 = np.where(finaldwarf.confidence_level == 1.0)[0]
conf2 = np.where(finaldwarf.confidence_level == 2.0)[0]
conf3 = np.where(finaldwarf.confidence_level == 3.0)[0]
print("Confidence Level -1 : {}".format(len(conf_1)), 
      "Confidence Level 0 : {}".format(len(conf0)), 
      "Confidence Level 1 : {}".format(len(conf1)), 
      'Confidence Level 2 : {}'.format(len(conf2)), 
      'Confidence Level 3 : {}'.format(len(conf3)))
if resolveflag:
    finalconf.to_csv('RESOLVE_snr5_master_conf_bary.csv')
if ecoflag:
    finalconf.to_csv('ECO_snr5_master_conf_bary.csv')
finalsfagn = finalconf.loc[unique]
#sfagn_targets = list(finalsfagn.iloc[np.where(finalsfagn['confidence_level'] > 0)].name)
sfagn_targets = [x for x in list(finalsfagn.name) if x[1] == 'f']
targetlist = master.loc[sfagn_targets][['radeg','dedeg']]
sami = ['rs0010']#, 'rs0775']
manga = ['rf0503']
#targetlist = targetlist.drop(manga)
#targetlist = targetlist.drop(sami)
'''


