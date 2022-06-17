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
#method = 'final'
#method = 'ap'
#method = 'exp'
#method = 'cog'


def create_df(data, datacols, dfcols):
    
    data_dict = {dfcols[0] : data[datacols[0]]}

    for i in range(1,len(dfcols)):
        data_dict[dfcols[i]] = data[datacols[i]]

    for x in data_dict.keys():
        if 'byteswap' in dir(data_dict[x]):
            data_dict[x] = data_dict[x].byteswap().newbyteorder() 

    df = pd.DataFrame(data = data_dict)
    df.index = df.name
    return df

def comparison_plot(df1, df2, band, df1matchcol, df2matchcol,
                    df1name, df2name, method, plotout):
    common, df1ndx, df2ndx = np.intersect1d(df1[df1matchcol], df2[df1matchcol],
                                            return_indices = True)
    if len(common) > 1:
        df1.index = df1[df1matchcol]
        df2.index = df2[df1matchcol]    
        residual_den = df1[band].iloc[df1ndx]
        residual_num = (df1[band].iloc[df1ndx]-df2[band].iloc[df2ndx])
    
    #   residual_num = (df2[band].loc[common])
        residual = residual_num#/residual_den
        residual_err_num = np.sqrt(df2['e'+band].iloc[df2ndx]**2 + 
                                   df1['e'+band].iloc[df1ndx]**2)
    #    residual_err_den = df1['e'+band].iloc[df1ndx]
        residual_err = residual_err_num#residual*(residual_err_num/residual_num + residual_err_den/residual_den)
        plt.figure()
        #plt.plot(reseco[band].loc[resecomatch],residual,'o')
        plt.errorbar(df1[band].iloc[df1ndx],residual,fmt='o',
                    yerr = np.array(residual_err))
    
    #    plt.errorbar(df1[band].loc[common],df2[band].loc[common],fmt='o',
    #                    xerr = np.array(df1['e'+band].loc[common]), 
    #                    yerr = np.array(df2['e'+band].loc[common]))
        
        plt.plot(np.arange(0,20), 0*np.arange(0,20))
        plt.plot(np.arange(0,20), 0*np.arange(0,20)+np.nanmedian(residual), 
                 'k--')
        plt.xlabel(df1name+' '+band)
        plt.ylabel('('+df1name+' '+band+' - '+df2name+' '+band+')')#/RESOLVE '+band)
        plt.xlim(min(df1[band].iloc[df1ndx])-0.25,
                 max(df1[band].iloc[df1ndx])+0.25)
        plt.title(method)
    #            plt.yscale('log')
    #            plt.xscale('log')
    #            plt.ylim(-2,2)
        if plotout:
            posndx = ((residual - residual_err)>0) & (residual>0 )  
            negndx = ((residual + residual_err)<0) & (residual<0 )  
            ndx = np.array(residual.index[posndx | negndx])
            plt.errorbar(df1[band].loc[ndx],residual.loc[ndx],fmt='o', color = 'red',
                        yerr = np.array(residual_err.loc[ndx]))
            return ndx
    else:
        return np.nan

wisecat = pd.read_csv('resolve_wise_db.csv')
wisecat.index = wisecat.name


methods = ['final']#, 'ap', 'exp', 'cog']
method = methods[0]

ecofile = readsav("../ecoSEL_wise_wisemask_031022.dat",python_dict=True)
ecodat = ecofile['resolve_wise'][method][0]
eco = create_df(ecodat, datacols = ['name', 'mw1', 'mw2', 'mw3', 'mw4', 'emw1', 'emw2', 'emw3', 'emw4'], 
                dfcols = ['name', 'mw1', 'mw2', 'mw3', 'mw4', 'emw1', 'emw2', 'emw3', 'emw4'])
eco = eco[(eco['mw1'] !=999) & (eco['mw1'] > 0)] 
ecofull = readsav("../eco_wresa_032918.dat")
eco['econame'] = np.zeros(len(eco))
ecomatch, ecofullndx, econdx = np.intersect1d(ecofull.names, eco.name, return_indices = True)
eco['econame'].iloc[econdx] = np.array(ecofull.econames[ecofullndx])
eco.index = eco.econame

resolve = readsav("../resolvecatalog_031622.dat",python_dict=True)
resphot = readsav("../resolvecatalogphot_031622.dat")

resfile = readsav("resolve_wise_wisemask_021522.dat",python_dict=True)
resdat = resfile['resolve_wise'][method][0]

res = create_df(resdat, datacols = ['name', 'mw1', 'mw2', 'mw3', 'mw4', 'emw1', 'emw2', 'emw3', 'emw4'], 
                dfcols = ['name', 'mw1', 'mw2', 'mw3', 'mw4', 'emw1', 'emw2', 'emw3', 'emw4'])
res['econame'] = resolve['econame']

#resphot['name'] = resolve['name']
#resphot['econame'] = resolve['econame']
#res = create_df(resphot, datacols = ['name', 'econame', 'mw1w', 'mw2w', 'mw3w', 'mw4w', 'emw1w', 'emw2w', 'emw3w', 'emw4w'], 
#                dfcols = ['name', 'econame', 'mw1', 'mw2', 'mw3', 'mw4', 'emw1', 'emw2', 'emw3', 'emw4'])
res = res[(res['mw1'] > 0) & (res['mw1'] !=999) ] 
res = res[np.isfinite(res['emw1'])] 

gamaecodf = pd.read_csv('GAMA_WISE_ECO_barysample.csv')
gamaecodf.index = gamaecodf.name
#Convert GAMA magnitudes to Vega to compare with general WISE mags
gamaecodf['mw1'] = gamaecodf['MAGPRO_W1'] - 2.699
gamaecodf['mw2'] = gamaecodf['MAGPRO_W2'] - 3.339
gamaecodf['mw3'] = gamaecodf['MAGPRO_W3'] - 5.174
gamaecodf['mw4'] = gamaecodf['MAGPRO_W4'] - 6.620
gamaecodf['emw1'] = gamaecodf['MAGPROERR_W1'] - 2.699
gamaecodf['emw2'] = gamaecodf['MAGPROERR_W2'] - 3.339
gamaecodf['emw3'] = gamaecodf['MAGPROERR_W3'] - 5.174
gamaecodf['emw4'] = gamaecodf['MAGPROERR_W4'] - 6.620
gamaecodf= gamaecodf[(gamaecodf['mw1'] !=999) & (gamaecodf['mw1'] > 0)] 

gamamethod = 'isophotal'
gamamethod = 'rec'
gamaresdf = pd.read_csv('GAMA_WISE_RESOLVE_barysample.csv')
gamaresdf.index = gamaresdf.name
#Convert GAMA magnitudes to Vega to compare with general WISE mags
if method == 'cog':
    gamaresdf['mw1'] = gamaresdf['MAGPRO_W1'] - 2.699
    gamaresdf['mw2'] = gamaresdf['MAGPRO_W2'] - 3.339
    gamaresdf['mw3'] = gamaresdf['MAGPRO_W3'] - 5.174
    gamaresdf['mw4'] = gamaresdf['MAGPRO_W4'] - 6.620
    gamaresdf['emw1'] = gamaresdf['MAGPROERR_W1'] - 2.699
    gamaresdf['emw2'] = gamaresdf['MAGPROERR_W2'] - 3.339
    gamaresdf['emw3'] = gamaresdf['MAGPROERR_W3'] - 5.174
    gamaresdf['emw4'] = gamaresdf['MAGPROERR_W4'] - 6.620
    gamaresdf= gamaresdf[(gamaresdf['mw1'] !=999) & (gamaresdf['mw1'] > 0)] 

elif method =='ap':
    #APERTURE PHOT
    gamaresdf['mw1'] = gamaresdf['MAGSTDAP_W1'] - 2.699
    gamaresdf['mw2'] = gamaresdf['MAGSTDAP_W2'] - 3.339
    gamaresdf['mw3'] = gamaresdf['MAGSTDAP_W3'] - 5.174
    gamaresdf['mw4'] = gamaresdf['MAGSTDAP_W4'] - 6.620
    gamaresdf['emw1'] = gamaresdf['MAGSTDAPERR_W1'] - 2.699
    gamaresdf['emw2'] = gamaresdf['MAGSTDAPERR_W2'] - 3.339
    gamaresdf['emw3'] = gamaresdf['MAGSTDAPERR_W3'] - 5.174
    gamaresdf['emw4'] = gamaresdf['MAGSTDAPERR_W4'] - 6.620

elif (method == 'final') & (gamamethod == 'isophotal'):
    gamaresdf['mw1'] = gamaresdf['MAGISO_W1'] - 2.699
    gamaresdf['mw2'] = gamaresdf['MAGISO_W2'] - 3.339
    gamaresdf['mw3'] = gamaresdf['MAGISO_W3'] - 5.174
    gamaresdf['mw4'] = gamaresdf['MAGISO_W4'] - 6.620
    gamaresdf['emw1'] = gamaresdf['MAGISOERR_W1'] - 2.699
    gamaresdf['emw2'] = gamaresdf['MAGISOERR_W2'] - 3.339
    gamaresdf['emw3'] = gamaresdf['MAGISOERR_W3'] - 5.174
    gamaresdf['emw4'] = gamaresdf['MAGISOERR_W4'] - 6.620
    gamaresdf= gamaresdf[(gamaresdf['mw1'] !=999) & (gamaresdf['mw1'] > 0)] 

elif method == 'final':
    gamaresdf['mw1'] = gamaresdf['MAG_W1'] - 2.699
    gamaresdf['mw2'] = gamaresdf['MAG_W2'] - 3.339
    gamaresdf['mw3'] = gamaresdf['MAG_W3'] - 5.174
    gamaresdf['mw4'] = gamaresdf['MAG_W4'] - 6.620
    gamaresdf['emw1'] = gamaresdf['MAGERR_W1'] - 2.699
    gamaresdf['emw2'] = gamaresdf['MAGERR_W2'] - 3.339
    gamaresdf['emw3'] = gamaresdf['MAGERR_W3'] - 5.174
    gamaresdf['emw4'] = gamaresdf['MAGERR_W4'] - 6.620
    gamaresdf= gamaresdf[(gamaresdf['mw1'] !=999) & (gamaresdf['mw1'] > 0)] 

reseco = res[res.econame != 'notineco']
reseco.index = reseco.econame
resecomatch, resndx, econdx = np.intersect1d(reseco['econame'], eco.econame, return_indices = True)
resecomatchcond = [x in resecomatch for x in eco.econame]



allwisecatorig = pd.read_csv('RESOLVE_WISE_allwise.csv')
#remove duplicates
unq = np.unique(allwisecatorig['name'])
allwisecat = pd.DataFrame({})
for i in range(len(unq)):
    df = allwisecatorig[allwisecatorig['name'] == unq[i]]
    if len(df):
        ndx = df.index.values[df['dist_x'] == min(df['dist_x'])]
        if len(ndx) !=1 :
            if allwisecatorig.name[ndx[0]] == allwisecatorig.name[ndx[1]]:
               ndx = ndx[0]
            else:
                ndx = []
            
        allwisecat = allwisecat.append(df.loc[ndx])


allwise = create_df(allwisecat, datacols = ['name', 'w1mpro', 'w2mpro', \
                                            'w3mpro', 'w4mpro', 'w1sigmpro', \
                                            'w2sigmpro', 'w3sigmpro', 'w4sigmpro'], 
                dfcols = ['name', 'mw1', 'mw2', 'mw3', 'mw4', 'emw1', 'emw2', 'emw3', 'emw4'])


wisecatorig = pd.read_csv('RESOLVE_WISE_allsky_orig.csv')
#remove duplicates
unq = np.unique(wisecatorig['name'])
wisecat = pd.DataFrame({})
for i in range(len(unq)):
    df = wisecatorig[wisecatorig['name'] == unq[i]]
    if len(df):
        ndx = df.index.values[df['dist_x'] == min(df['dist_x'])]
        if len(ndx) !=1 :
            if wisecatorig.name[ndx[0]] == wisecatorig.name[ndx[1]]:
               ndx = ndx[0]
            else:
                ndx = []
            
        wisecat = wisecat.append(df.loc[ndx])


wiseallsky = create_df(wisecat, datacols = ['name', 'w1mpro', 'w2mpro', \
                                            'w3mpro', 'w4mpro', 'w1sigmpro', \
                                            'w2sigmpro', 'w3sigmpro', 'w4sigmpro'], 
                dfcols = ['name', 'mw1', 'mw2', 'mw3', 'mw4', 'emw1', 'emw2', 'emw3', 'emw4'])


bands = ['mw1']#, 'mw2','mw3','mw4']
for band in bands:
#    outliers = comparison_plot(res, eco, band, 'econame', 'econame',
#                    'RESOLVE DB', 'NEW ECO', method = 'final', plotout = 1)        

#    outliers = comparison_plot(gamaecodf, eco, band, df1matchcol = 'name', 
#                               df2matchcol = 'econame', df1name = 'GAMA WISE', 
#                               df2name = 'NEW ECO', method = 'final', plotout = 1)        

#    outliers = comparison_plot(gamaresdf, res, band, df1matchcol = 'name', 
#                               df2matchcol = 'name', df1name = 'GAMA isophotal', 
#                               df2name = 'RESOLVE DB', method = method, plotout = 1)        
    
#    outliers = comparison_plot(allwise, res, band, df1matchcol = 'name', 
#                               df2matchcol = 'name', df1name = 'ALLWISE profile-fit', 
#                               df2name = 'RESOLVE DB', method = method, plotout = 1)        
    
#    outliers = comparison_plot(wiseallsky, res, band, df1matchcol = 'name', 
#                               df2matchcol = 'name', df1name = 'WISE All-Sky profile-fit', 
#                               df2name = 'RESOLVE DB', method = method, plotout = 1)        
    
#    outliers = comparison_plot(gamaresdf, allwise, band, df1matchcol = 'name', 
#                               df2matchcol = 'name', df1name = 'GAMA recommended', 
#                               df2name = 'ALLWISE profile-fit', method = method, plotout = 1)        
    
    outliers = comparison_plot(gamaresdf, wiseallsky, band, df1matchcol = 'name', 
                               df2matchcol = 'name', df1name = 'GAMA recommended', 
                               df2name = 'WISE-allsky profile-fit', method = method, plotout = 1)        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
#df = pd.read_csv("../sfr_nuv_wisew.txt")
##                 names = ['name','mw3w','mw4w','mw1w','groupcz','deextrestnuvmag','sfr_nuv','sfr_nuv_wisew'])
#df.index = df.econame
#df = df[df.mw3w != -999]

    
        #bands = ['mw1']#, 'mw2','mw3','mw4']
        #for band in bands:
        #    
        ##    plt.figure()
        ##    plt.title(band)
        ##    plt.errorbar(reseco[band].loc[resecomatch],eco[band].loc[resecomatch],fmt='o',
        ##                xerr = np.array(reseco['e'+band].loc[resecomatch]), 
        ##                yerr = np.array(eco['e'+band].loc[resecomatch]))
        ##    plt.plot(np.arange(0,20), np.arange(0,20))
        ##    plt.xlabel('RESOLVE mags')
        ##    plt.ylabel('ECO mags')
        ##    plt.xlim(5,17.5)
        ##    plt.ylim(5,17.5)
        #    
        #    residual_den = reseco[band].loc[resecomatch]
        #    residual_num = (reseco[band].loc[resecomatch]-eco[band].loc[resecomatch])
        #    residual = residual_num#/residual_den
        #    residual_err_num = np.sqrt(reseco['e'+band].loc[resecomatch]**2 + eco['e'+band].loc[resecomatch]**2)
        #    residual_err_den = reseco['e'+band].loc[resecomatch]
        #    residual_err = residual_err_num#residual*(residual_err_num/residual_num + residual_err_den/residual_den)
        #    plt.figure()
        #    #plt.plot(reseco[band].loc[resecomatch],residual,'o')
        #    plt.errorbar(reseco[band].loc[resecomatch],residual,fmt='o',
        #                yerr = np.array(residual_err))
        #    #            xerr = np.array(reseco['e'+band].loc[resecomatch]), 
        #    
        #    plt.plot(np.arange(0,20), 0*np.arange(0,20))
        #    plt.plot(np.arange(0,20), 0*np.arange(0,20)+np.nanmedian(residual), 'k--')
        #    plt.xlabel('RESOLVE '+band)
        #    plt.ylabel('(RESOLVE '+band+' - ECO '+band+')')#/RESOLVE '+band)
        #    plt.xlim(min(reseco[band].loc[resecomatch])-0.25,max(reseco[band].loc[resecomatch])+0.25)
        #    plt.ylim(-2,2)
        #    posndx = ((residual - residual_err)>0) & (residual>0 )  
        #    negndx = ((residual + residual_err)<0) & (residual<0 )  
        #    ndx = np.array(residual.index[posndx | negndx])
        #    plt.errorbar(reseco[band].loc[ndx],residual.loc[ndx],fmt='o', color = 'red',
        #                yerr = np.array(residual_err.loc[ndx]))
        

#resbary = pd.read_csv('../RESOLVE_barysample.csv')
#for i in range(np.where(resbary.name == 'rf0336')[0]):#len(resbary)):
#    if resbary.dedeg[i] > 0:
#        txt = "| {ra:3.8f} | {dec:2.8f} | {name} | "
##        txt = "  {ra:3.8f}   {dec:2.8f}  "
#        print(txt.format(ra= resbary.radeg[i], dec = resbary.dedeg[i], 
#                         name = resbary.name[i]))
#    else:
#        txt = "| {ra:3.8f} | {dec:2.7f} | {name} | "
##        txt = "  {ra:3.8f}   {dec:2.7f}  "
#        print(txt.format(ra= resbary.radeg[i], dec = resbary.dedeg[i], 
#                         name = resbary.name[i]))