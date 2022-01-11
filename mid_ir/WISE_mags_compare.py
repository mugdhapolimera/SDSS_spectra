# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 16:19:43 2020

@author: mugdhapolimera
"""

import numpy as np
from scipy.io.idl import readsav
import pandas as pd
import matplotlib.pyplot as plt
#resfile = readsav("resolve_wise_102919.dat") 
method = 'ap'#'exp'#'cog'

resolve = readsav("resolvecatalog.dat",python_dict=True)
resphot = readsav("resolvecatalogphot.dat")

ecofile = readsav("ecowresa_wise_wisemask_resvalues_111621.dat",python_dict=True)
ecodat = ecofile['resolve_wise'][method][0]
eco = {'name' : ecodat['name'],
               'mw1' : ecodat['mw1'],
              'mw2' : ecodat['mw2'],
              'mw3' : ecodat['mw3'],
              'mw4' : ecodat['mw4'],
              'emw1' : ecodat['emw1'],
              'emw2' : ecodat['emw2'],
              'emw3' : ecodat['emw3'],
              'emw4' : ecodat['emw4']}

for x in eco.keys():
    eco[x] = eco[x].byteswap().newbyteorder() 

eco = pd.DataFrame(data = eco)
eco.index = eco.name

eco = eco[(eco['mw1'] !=999) & (eco['mw1'] > 0)] 
#eco = eco[eco['mw1']/eco['emw1'] > 5]

ecofull = readsav("../eco_wresa_032918.dat")

eco['econame'] = np.zeros(len(eco))

ecomatch, ecofullndx, econdx = np.intersect1d(ecofull.names, eco.name, return_indices = True)
eco['econame'].iloc[econdx] = ecofull.econames[ecofullndx]
eco['nuvmag'] = 1
eco['nuvmag'].iloc[econdx] = ecofull.nuvmag[ecofullndx]

resolve_wisemask_file = readsav("resolve_wise_wisemask_010320.dat", python_dict = True)
resolve_wisemask = resolve_wisemask_file['resolve_wise'][method][0]

#res = {'name' : resolve['name'],
#       'econame' : resolve['econame'],
#              'mw1' : resphot['mw1w'],
#              'mw2' : resphot['mw2w'],
#              'mw3' : resphot['mw3w'],
#              'mw4' : resphot['mw4w'],
#              'emw1' : resphot['emw1w'],
#              'emw2' : resphot['emw2w'],
#              'emw3' : resphot['emw3w'],
#              'emw4' : resphot['emw4w']}
#res = {'name' : resolve['name'],
#       'econame' : resolve['econame'],
#              'mw1' : resphot['mw1o'],
#              'mw2' : resphot['mw2o'],
#              'mw3' : resphot['mw3o'],
#              'mw4' : resphot['mw4o'],
#              'emw1' : resphot['emw1o'],
#              'emw2' : resphot['emw2o'],
#              'emw3' : resphot['emw3o'],
#              'emw4' : resphot['emw4o']}
res = {'name' : resolve['name'],
       'econame' : resolve['econame'],
       'sfr_nuv_wise': resolve['sfr_nuv_wise'],
       'sfr_nuv': resolve['sfr_nuv'],
       'nuvmag': resphot['nuvmag'],
               'mw1' : resolve_wisemask['mw1'],
              'mw2' : resolve_wisemask['mw2'],
              'mw3' : resolve_wisemask['mw3'],
              'mw4' : resolve_wisemask['mw4'],
              'emw1' : resolve_wisemask['emw1'],
              'emw2' : resolve_wisemask['emw2'],
              'emw3' : resolve_wisemask['emw3'],
              'emw4' : resolve_wisemask['emw4']}

#Changing byter order from IDL dtype format to pandas default dtype
for x in res.keys():
    res[x] = res[x].byteswap().newbyteorder() 

#Making a pandas dataframe using the dictionary
res = pd.DataFrame(data = res)
res.index = res.name
res = res[res['mw1'] !=0] 
res = res[np.isfinite(res['emw1'])] 
reseco = res[res.econame != 'notineco']
reseco.index = reseco.econame
resecomatch, resndx, econdx = np.intersect1d(reseco['econame'], eco.econame, return_indices = True)
resecomatchcond = [x in resecomatch for x in eco.econame]
eco.index = eco.econame

bands = ['mw1']#, 'mw2','mw3','mw4']
for band in bands:
    
#    plt.figure()
#    plt.title(band)
#    plt.errorbar(reseco[band].loc[resecomatch],eco[band].loc[resecomatch],fmt='o',
#                xerr = np.array(reseco['e'+band].loc[resecomatch]), 
#                yerr = np.array(eco['e'+band].loc[resecomatch]))
#    plt.plot(np.arange(0,20), np.arange(0,20))
#    plt.xlabel('RESOLVE mags')
#    plt.ylabel('ECO mags')
#    plt.xlim(5,17.5)
#    plt.ylim(5,17.5)
    
    residual_den = reseco[band].loc[resecomatch]
    residual_num = (reseco[band].loc[resecomatch]-eco[band].loc[resecomatch])
    residual = residual_num#/residual_den
    residual_err_num = np.sqrt(reseco['e'+band].loc[resecomatch]**2 + eco['e'+band].loc[resecomatch]**2)
    residual_err_den = reseco['e'+band].loc[resecomatch]
    residual_err = residual_err_num#residual*(residual_err_num/residual_num + residual_err_den/residual_den)
    plt.figure()
    #plt.plot(reseco[band].loc[resecomatch],residual,'o')
    plt.errorbar(reseco[band].loc[resecomatch],residual,fmt='o',
                yerr = np.array(residual_err))
    #            xerr = np.array(reseco['e'+band].loc[resecomatch]), 
    
    plt.plot(np.arange(0,20), 0*np.arange(0,20))
    plt.plot(np.arange(0,20), 0*np.arange(0,20)+np.nanmedian(residual), 'k--')
    plt.xlabel('RESOLVE '+band)
    plt.ylabel('(RESOLVE '+band+' - ECO '+band+')')#/RESOLVE '+band)
    plt.xlim(min(reseco[band].loc[resecomatch])-0.25,max(reseco[band].loc[resecomatch])+0.25)
    plt.ylim(-2,2)
    posndx = ((residual - residual_err)>0) & (residual>0 )  
    negndx = ((residual + residual_err)<0) & (residual<0 )  
    ndx = np.array(residual.index[posndx | negndx])
    plt.errorbar(reseco[band].loc[ndx],residual.loc[ndx],fmt='o', color = 'red',
                yerr = np.array(residual_err.loc[ndx]))


df = pd.read_csv("../sfr_nuv_wisew.txt")
#                 names = ['name','mw3w','mw4w','mw1w','groupcz','deextrestnuvmag','sfr_nuv','sfr_nuv_wisew'])
df.index = df.econame
df = df[df.mw3w != -999]
resdfmatch, resndx, dfndx = np.intersect1d(reseco['econame'], df.econame, return_indices = True)
resdfmatchcond = [x in resdfmatch for x in df.econame]

residual_den = reseco['sfr_nuv_wise'].loc[resecomatch]
residual_num = (reseco['sfr_nuv_wise'].loc[resecomatch]-df['sfr_nuv_wisew'].loc[resecomatch])
residual = residual_num/residual_den
plt.figure()
plt.plot(reseco['sfr_nuv_wise'].loc[resecomatch],residual,'o')
#plt.plot(np.arange(0,20), 0*np.arange(0,20))
#plt.plot(np.arange(0,20), 0*np.arange(0,20)+np.nanmedian(residual), 'k--')
plt.xlabel('RESOLVE SFR_NUV_WISE')
plt.ylabel('SFR_NUV_WISE residual % (RESOLVE -ECO)/RESOLVE)')
#plt.yscale("log")
#plt.xlim(min(reseco[band].loc[resecomatch])-0.25,max(reseco[band].loc[resecomatch])+0.25)
#plt.ylim(-2,2)
plt.figure()
plt.plot(reseco['mw1'].loc[resecomatch],residual,'o')
plt.xlabel('RESOLVE MW1')
plt.ylabel('SFR_NUV_WISE residual % (RESOLVE -ECO)/RESOLVE)')
#plt.yscale("log")
plt.xlim(16,7)

residual_den = reseco['sfr_nuv'].loc[resecomatch]
residual_num = (reseco['sfr_nuv'].loc[resecomatch]-df['sfr_nuv'].loc[resecomatch])
residual = residual_num#/residual_den
plt.figure()
plt.plot(reseco['sfr_nuv'].loc[resecomatch],residual,'o')
#plt.plot(np.arange(0,20), 0*np.arange(0,20))
#plt.plot(np.arange(0,20), 0*np.arange(0,20)+np.nanmedian(residual), 'k--')
plt.xlabel('RESOLVE SFR_NUV')
plt.ylabel('(RESOLVE SFR_NUV - ECO SFR_NUV)')#/RESOLVE '+band)

#resecoinobs = [x for x in np.array(resecomatch) if resolve['inobssample'][np.where(resolve['name'] == reseco.name.loc[x])]]
#plt.plot(reseco['mw1'].loc[resecoinobs],residual.loc[resecoinobs],'ro', mfc = 'none')


residual_den = reseco['nuvmag'].loc[resecomatch]
residual_num = (reseco['nuvmag'].loc[resecomatch]-eco['nuvmag'].iloc[econdx])
residual = residual_num#/residual_den
#residual = residual[residual != 1]
plt.figure()
plt.plot(reseco['nuvmag'].loc[resecomatch],residual,'o')
#plt.plot(np.arange(0,20), 0*np.arange(0,20))
#plt.plot(np.arange(0,20), 0*np.arange(0,20)+np.nanmedian(residual), 'k--')
plt.xlabel('RESOLVE NUV mag')
plt.ylabel('NUV mag residual (RESOLVE -ECO)')
plt.xlim(23,14)

resdf = pd.read_csv('../RESOLVElive2022Jan7.csv')
resdf.index = resdf.name
ecodf = pd.read_csv('../ECOlive2022Jan7.csv')
ecodf.index = ecodf.name
resecodf = resdf[resdf.econame != 'notineco']
resecodf.index = resecodf.econame
resecodfmatch, resdfndx, ecodfndx = np.intersect1d(resecodf['econame'], 
                                        ecodf.name, return_indices = True)

residual_den = resecodf['nuvmag'].loc[resecodfmatch]
residual_num = (resecodf['nuvmag'].loc[resecodfmatch] - ecodf['extnuv'].loc[resecodfmatch])
residual = residual_num#/residual_den
#residual = residual[residual != 1]
plt.figure()
plt.plot(resecodf['nuvmag'].loc[resecodfmatch],residual,'o')
#plt.plot(np.arange(0,20), 0*np.arange(0,20))
#plt.plot(np.arange(0,20), 0*np.arange(0,20)+np.nanmedian(residual), 'k--')
plt.xlabel('RESOLVE NUV mag')
plt.ylabel('NUV mag residual (RESOLVE -ECO)')
plt.xlim(23,14)

residual_den = resecodf['sfr_nuv'].loc[resecodfmatch]
residual_num = (resecodf['sfr_nuv'].loc[resecodfmatch] - ecodf['sfr_nuv'].loc[resecodfmatch])
residual = residual_num#/residual_den
#residual = residual[residual != 1]
plt.figure()
plt.plot(resecodf['nuvmag'].loc[resecodfmatch],residual,'o')
#plt.plot(np.arange(0,20), 0*np.arange(0,20))
#plt.plot(np.arange(0,20), 0*np.arange(0,20)+np.nanmedian(residual), 'k--')
plt.xlabel('RESOLVE NUV mag')
plt.ylabel('NUV mag residual (RESOLVE -ECO)')
plt.xlim(23,14)









residual_den = resecodf['nuvmag'].loc[resecomatch]
residual_num = (resecodf['nuvmag'].loc[resecomatch] - eco['nuvmag'].iloc[econdx])
residual = residual_num#/residual_den
#residual = residual[residual != 1]
plt.figure()
plt.plot(resecodf['nuvmag'].loc[resecomatch],residual,'o')
#plt.plot(np.arange(0,20), 0*np.arange(0,20))
#plt.plot(np.arange(0,20), 0*np.arange(0,20)+np.nanmedian(residual), 'k--')
plt.xlabel('RESOLVE NUV mag')
plt.ylabel('NUV mag residual (RESOLVE -ECO)')
plt.xlim(23,14)


