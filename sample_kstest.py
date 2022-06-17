# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 16:25:48 2021

@author: mugdhapolimera
"""
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import matplotlib
from scipy import stats
matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams.update({'axes.linewidth': 2})
matplotlib.rcParams.update({'lines.linewidth': 2})

res = 0
eco = 1
if res:
    inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_blend_dext_new.csv'
    df = pd.read_csv(inputfile)
    df.index = df.name

    print( len(df))
    ra=df.radeg
    dec=df.dedeg
    grpcz = df.grpcz
    cz = df.cz
    infall = (ra > 22*15.) | (ra < 3*15.)
    inspring = (ra > 8.75*15.) & (ra < 15.75*15.)
    inobssample = ((grpcz >= 4500.) & (grpcz <= 7000.) & (infall | inspring))# | inspring)) 
    df = df[inobssample]
    parent = df
    flinsample = df.fl_insample
    mgas = df.logmgas
    mstars = df.logmstar
    mbary = 10**mgas + 10**mstars
    resfullmbary = np.log10(mbary)
    inobssample = (resfullmbary > 9.2)
    res_el_full = pd.read_csv('RESOLVE_full_hasnr5_dext_jhu.csv')
    inspring = (res_el_full.radeg > 8.75*15.) & (res_el_full.radeg < 15.75*15.)
    infall = (res_el_full.radeg > 22*15.) | (res_el_full.radeg < 3*15.)
#    res_el = res_el[inspring]
#    res_el = res_el[infall]
    res_sel = pd.read_csv('RESOLVE_full_snr5_dext_jhu.csv')
    
    inspring = (res_sel.radeg > 8.75*15.) & (res_sel.radeg < 15.75*15.)
    infall = (res_sel.radeg > 22*15.) | (res_sel.radeg < 3*15.)
#    res_sel = res_sel[inspring]
#    res_sel = res_sel[infall]
    jhu = pd.read_csv('RESOLVE_full_snr5_dext_jhu.csv')
#    jhu = jhu[(jhu.radeg > 8.75*15.) & (jhu.radeg < 15.75*15.)]
#    jhu = jhu[(jhu.radeg > 22*15.) | (jhu.radeg < 3*15.)]

    port = pd.read_csv('RESOLVE_full_snr5_dext_port.csv')
#    port = port[(port.radeg > 8.75*15.) & (port.radeg < 15.75*15.)]
#    port = port[(port.radeg > 22*15.) | (port.radeg < 3*15.)]

    nsa = pd.read_csv('RESOLVE_full_snr5_dext_nsa.csv')
#    nsa = nsa[(nsa.radeg > 8.75*15.) & (nsa.radeg < 15.75*15.)]
#    nsa = nsa[(nsa.radeg > 22*15.) | (nsa.radeg < 3*15.)]

if eco:
    inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_full_blend_dext_new.csv'
    df = pd.read_csv(inputfile)
    df.index = df.name
    print (len(df))
    mgas = df.logmgas
    mstars = df.logmstar
    mbary = 10**mgas + 10**mstars
    ineco = (130.05 < df.radeg) & (df.radeg < 237.45)
    inobssample = (((df.grpcz >= 3000.) & (df.grpcz <= 7000.)) & 
                   (np.log10(mbary) > 9.2) & ineco)#\
    res_el_full = pd.read_csv('ECO_full_hasnr5_dext_jhu.csv')
    res_sel = pd.read_csv('ECO/SEL/ECO_full_snr5_dext_jhu.csv')
    res_el = res_el_full.loc[~res_el_full.index.isin(list(res_sel.index))]
   
    jhu = pd.read_csv('ECO/SEL/ECO_full_snr5_dext_jhu.csv')
    port = pd.read_csv('ECO/SEL/ECO_full_snr5_dext_port.csv')
    nsa = pd.read_csv('ECO/SEL/ECO_full_snr5_dext_nsa.csv')

res_sel.index = res_sel.name
res_el_full.index = res_el_full.name
    
resolve = df[inobssample]
resolve['logmbary'] = np.log10(10**resolve.logmstar + 10**resolve.logmgas)

res_el_full['logmbary'] = np.log10(10**res_el_full.logmstar + 10**res_el_full.logmgas)
res_sel['logmbary'] = np.log10(10**res_sel.logmstar + 10**res_sel.logmgas)

res_el = res_el_full.loc[~res_el_full.index.isin(list(res_sel.index))]
nonel_resolve = resolve.loc[~resolve.index.isin(list(res_el_full.index))]
nonsel_resolve = resolve.loc[~resolve.index.isin(list(res_sel.index))]


props = ['logmbary', 'logmh', 'modelu_rcorr']

for prop in props:

#    DD, pnullks1 = stats.ks_2samp(nonel_resolve[prop],res_el_full[prop])
    DD, pnullks1 = stats.ks_2samp(resolve[prop],res_el_full[prop], mode = 'asymp')
    print('Mass-limited & full-EL samples '+prop+' K-S test p-value = ' + str(pnullks1))
    
#    DD, pnullks2 = stats.ks_2samp(nonsel_resolve[prop], res_sel[prop])
    DD, pnullks2 = stats.ks_2samp(resolve[prop], res_sel[prop], mode = 'asymp')
    print('Mass-limited & SEL samples '+prop+' K-S test p-value = ' + str(pnullks2))
    
#    DD, pnullks3 = stats.ks_2samp(res_el[prop], res_sel[prop])
    DD, pnullks3 = stats.ks_2samp(res_el_full[prop], res_sel[prop], mode = 'asymp')
    print('SEL & non-sel EL samples '+prop+' K-S test p-value = ' + str(pnullks3)+'\n\n')
    



#bins = np.arange(0.5,3.5,0.25)
#plt.figure()
#plt.hist(resolve['modelu_rcorr'], bins = bins, color = 'gray', histtype = 'step', lw = 3, density = True)    
#plt.hist(res_el['modelu_rcorr'], bins = bins, color = 'orange', histtype = 'step', lw = 3, density = True)    
#plt.hist(res_sel['modelu_rcorr'], bins = bins, color = 'black', histtype = 'step', lw = 3, density = True)    
#
#bins = np.arange(10,14,0.5)
#plt.figure()
#plt.hist(resolve['logmh'], bins = bins, color = 'gray', histtype = 'step', lw = 3, density = True)    
#plt.hist(res_el['logmh'], bins = bins, color = 'orange', histtype = 'step', lw = 3, density = True)    
#plt.hist(res_sel['logmh'], bins = bins, color = 'black', histtype = 'step', lw = 3, density = True)    
#
#bins = np.arange(9.2,11,0.2)
#plt.figure()
#plt.hist(resolve['logmbary'], bins = bins, color = 'gray', histtype = 'step', lw = 3, density = True)    
#plt.hist(res_el['logmbary'], bins = bins, color = 'orange', histtype = 'step', lw = 3, density = True)    
#plt.hist(res_sel['logmbary'], bins = bins, color = 'black', histtype = 'step', lw = 3, density = True)    
#

jhu['logmbary'] = np.log10(10**jhu.logmstar + 10**jhu.logmgas)
port['logmbary'] = np.log10(10**port.logmstar + 10**port.logmgas)
nsa['logmbary'] = np.log10(10**nsa.logmstar + 10**nsa.logmgas)

bins = {'logmbary': np.arange(9.2,11,0.2), 
        'logmh': np.arange(10,14,0.5), 
        'modelu_rcorr': np.arange(0.5,3.5,0.25)}
# for prop in props:
#     print(prop)
# #    DD, pnullks = stats.ks_2samp(jhu[prop],port[prop], mode = 'asymp')
# #     DD, pnullks = stats.mannwhitneyu(jhu[prop],port[prop], alternative = 'two-sided')
# #    statistic,crit, pnullks = stats.anderson_ksamp([port[prop], jhu[prop]])
# #    confidence=stats.norm.interval(1-pnullks) # you fill in ???
#     DD, pnullks = stats.epps_singleton_2samp(jhu[prop],port[prop])
#     print('JHU and Port samples  K-S test p-value = ' + str(round(pnullks,4)))#[0:4])
# #    print('JHU and Port samples  K-S test confidence level = ' + str(confidence[1]))#[0:4])
    
# #    DD, pnullks = stats.ks_2samp(nsa[prop],port[prop], mode = 'asymp')
#     # DD, pnullks = stats.mannwhitneyu(port[prop], nsa[prop], alternative = 'two-sided')
# #    statistic,crit, pnullks = stats.anderson_ksamp([nsa[prop], port[prop]])
# #    confidence=stats.norm.interval(1-pnullks) # you fill in ???
#     DD, pnullks = stats.epps_singleton_2samp(nsa[prop],port[prop])
#     print('Port and NSA samples  K-S test p-value = ' + str(round(pnullks,4)))#[0:4])
# #    print('Port and NSA samples  K-S test confidence level = ' + str(confidence[1]))#[0:4])
    
# #    DD, pnullks = stats.ks_2samp(jhu[prop],nsa[prop], mode = 'asymp')
#     # DD, pnullks = stats.mannwhitneyu(nsa[prop], jhu[prop], alternative = 'two-sided')
# #    statistic,crit, pnullks = stats.anderson_ksamp([nsa[prop], jhu[prop]]) #significance level is stored in pnullks for ease of printing
# #    confidence=stats.norm.interval(1-pnullks) # you fill in ???
#     DD, pnullks = stats.epps_singleton_2samp(jhu[prop],nsa[prop])
#     print('NSA and JHU samples  K-S test p-value = ' + str(round(pnullks,4)))
# #    print('NSA and JHU samples  K-S test confidence level = ' + str(confidence[1])+'\n\n')#[0:4])
    
#     density = True
#     plt.figure()
#     plt.hist(jhu[prop], bins = bins[prop], color = 'blue', histtype = 'step', 
#              lw = 3, density = density, label = 'JHU SEL')    
#     plt.hist(port[prop], bins = bins[prop], color = 'orange', histtype = 'step', 
#              hatch = '//', lw = 3, density = density, label = 'Portsmouth SEL')    
#     plt.hist(nsa[prop], bins = bins[prop], color = 'black', histtype = 'step', 
#              hatch = '\\', lw = 3, density = density, label = 'NSA SEL')  
#     plt.legend()
#     plt.xlabel(prop)
#     plt.ylabel('Normalized # of galaxies')  

# plt.figure()
# plt.plot(jhu.radeg, jhu.dedeg, 'o', color = 'blue', mew = 0, alpha = 0.3, label = 'JHU SEL')    
# plt.plot(port.radeg, port.dedeg, 'o', color = 'orange', mew = 0, alpha = 0.3, label = 'JHU SEL')    
# plt.plot(nsa.radeg, nsa.dedeg, 'o', color = 'black', mew = 0, alpha = 0.3, label = 'JHU SEL')    
# plt.legend()





#import scipy.stats as st
#samp_a = np.sort(resmbary_el)
#samp_b = np.sort(resmbary_sel)
#samp_conc = np.sort(np.concatenate((samp_a, samp_b)))
#samp_a_cdf = [np.round(st.percentileofscore(samp_a, value)/100, 1) for value in samp_conc]
#
##cdf of sample b
#samp_b_cdf = [np.round(st.percentileofscore(samp_b, value)/100, 1) for value in samp_conc]
##print('samp_b_cdf:\n{}'.format(samp_b_cdf))
##print(20*'-')
#
##compute absolute difference
#samp_diff = np.abs(np.subtract(samp_a_cdf, samp_b_cdf))
#print('samp_diff:\n{}'.format(samp_diff))
#plt.figure(); 
#plt.plot(samp_conc,samp_a_cdf)
#plt.plot(samp_conc,samp_b_cdf)
#plt.plot(samp_conc,min(samp_a_cdf,samp_b_cdf)+samp_diff)

#print("D_crit ", 1.36*np.sqrt(1.0/len(samp_a) + 1.0/len(samp_b)))
#print("D_n ", max(samp_diff))
    
bins = np.arange(9.2,11,0.2)
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)
#plt.hist(np.log10(mbary), density = True, bins = bins, color = 'gainsboro', 
#         label = 'Parent Sample')
ax1.hist(nonel_resolve.logmbary, density = True, bins = bins,\
         color = 'gainsboro', lw = 2, label = 'Mass-limited Parent Sample')
ax1.hist(res_el.logmbary, density = True, histtype = 'step', bins = bins,
         color = 'orange', lw = 4, label = 'Emission Line Sample')
ax1.hist(res_sel.logmbary, density = True, histtype = 'step', bins = bins,
         color = 'black', lw = 3, hatch = 'X', label = 'Strong Emission Line Sample')
#plt.axvline(9.2)
ax1.set_xlim(9.0,11.2)
ax1.set_xlabel(r"$\rm log(M_{baryonic}/M_\odot)$")
ax1.set_ylabel("Normalized # of galaxies")
#plt.legend(loc = "upper right", fontsize = 15)
ax1.set_xticks(bins)
ax1.tick_params(labelsize = 14)
ax2.tick_params(labelsize = 14)
ax3.tick_params(labelsize = 14)

bins = np.arange(10,14,0.25)
#plt.hist(parent.logmh, density = True, bins = bins, color = 'gainsboro', 
#         label = 'Parent Sample')
ax2.hist(nonel_resolve.logmh, density = True, bins = bins,\
         color = 'gainsboro', lw = 2, label = 'Mass-limited Parent Sample')
ax2.hist(res_el.logmh, density = True, histtype = 'step', bins = bins,
         color = 'orange', lw = 4, label = 'Emission Line Sample')
ax2.hist(res_sel.logmh, density = True, histtype = 'step', bins = bins,
         color = 'black', lw = 3, hatch = 'X', label = 'Strong Emission Line Sample')
#plt.axvline(9.2)
ax2.set_xlim(9.75,14.)
ax2.set_xticks(np.arange(10,14,0.5))
ax2.set_xlabel(r"$\rm log(M_{halo}/M_\odot)$")
#ax2.set_ylabel("Normalized # of galaxies")
#ax2.legend(loc = "upper right", fontsize = 15)
ax1.tick_params(length=8, width=2)
ax2.tick_params(length=8, width=2)

bins = np.arange(0.5,3.5,0.25)
ax3.hist(nonel_resolve.modelu_rcorr, density = True, bins = bins,\
         color = 'gainsboro', lw = 2, label = 'Mass-limited Parent Sample')
ax3.hist(res_el.modelu_rcorr, density = True, histtype = 'step', bins = bins,
         color = 'orange', lw = 4, label = 'Emission Line Sample')
ax3.hist(res_sel.modelu_rcorr, density = True, histtype = 'step', bins = bins,
         color = 'black', lw = 3, hatch = 'X', label = 'Strong Emission Line Sample')
#plt.axvline(9.2)
ax3.set_xlim(0.25,3.5)
ax3.set_xticks(bins)
ax3.set_xlabel(r"$\rm (u-r)^e$ color")
#ax3.set_ylabel("Normalized # of galaxies")
ax3.legend(loc = "upper right", fontsize = 15)
ax3.tick_params(length=8, width=2)