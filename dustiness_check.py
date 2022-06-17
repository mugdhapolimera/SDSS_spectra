# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 14:07:04 2021

@author: mugdhapolimera
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("ECO+RESOLVE_snr5_jhu.csv")#pd.read_csv("ECO+RESOLVE_snr5_jhu.csvEBV_excess_.csv")
df.index = df.name
df = df.drop(index = 'ECO03033')
df['EBV_excess'] = df['h_alpha_flux']/df['h_beta_flux']/2.86

jhuflag = pd.read_csv('ECO/SEL/eco+resolve_emlineclass_dext_snr5_jhu.csv')
jhuflag.index = jhuflag.galname
jhuflag = jhuflag.drop(index = 'ECO03033')

jhu = pd.read_csv('ECO+RESOLVE_snr5_dext_jhu.csv')#[jhuflag.sftoagn]
jhu.index = jhu.name
jhu = jhu.drop(index = 'ECO03033')
print np.nanmedian(df.EBV_excess), np.nanmedian(df.EBV_excess[jhu.logmstar < 9.5]),\
 np.nanmedian(df.EBV_excess[jhu.logmstar > 9.5])

plt.figure()
#plt.plot(jhu.logmstar, df.EBV_excess,'.', color = 'gray')
#plt.plot(jhu.logmstar[jhuflag.agntosf], df.EBV_excess[jhuflag.agntosf],'c^', ms = 15)
#for i, txt in enumerate(jhuflag.galname[jhuflag.agntosf]):
#    plt.annotate(txt, (jhu.logmstar.loc[jhuflag.galname[jhuflag.agntosf][i]], 
#                             df.EBV_excess.loc[jhuflag.galname[jhuflag.agntosf][i]]))
keys = ['defstarform', 'defagn', 'composite', 'sftoagn','agntosf']
    
marker = {'agntosf': 'c^', 'ambigagn': 'ms', 'composite': 'ms', 'defagn': 'rs', 
          'defliner': 'yo', 'defseyf': 'co', 'heiisel': 'ks',
          'defstarform': 'k.','sftoagn': 'bs', 'sftoagn1': 's', 'sftoagn2': 'm*'}

colors = {'agntosf': 'c', 'ambigagn': 'm', 'composite': 'm', 'defagn': 'r', 
          'defliner': 'y', 'defseyf': 'c', 'heiisel': 'k',
          'defstarform': 'gray', 'sftoagn': 'b', 'sftoagn2': 'b'}

labels = {'agntosf': 'Low-SII AGN', 'ambigagn': 'Ambiguous AGN', 
          'composite': 'Composite', 'defagn': 'Definite AGN', 
          'defliner': 'LINER', 'defseyf': 'Seyfert', 
          'heiisel': 'HeII-Selected AGN', 'defstarform': 'Definite SF', 
          'sftoagn': 'SFing-AGN', 'sftoagn2' : 'MP-AGN2'}

percent = {'agntosf': 0, 'composite': 0, 'defagn': 0, 'heiisel': 57,
          'defstarform': 0, 'sftoagn': 0}

for key in keys:
    sel = np.array(df.loc[jhuflag[key]].name)
    if key == 'agntosf':
        zo = 6
        plt.plot(df.logmstar.loc[sel], df.EBV_excess.loc[sel],
            marker[key], ms = 14, label = labels[key], zorder = zo, mec = 'k', mew = 2)

    else:
        zo = 2
    if key == 'defstarform':
        plt.plot(df.logmstar.loc[sel], df.EBV_excess.loc[sel],
            '.', color = 'gray', label = labels[key], zorder = zo, alpha = 0.5, mew = 0)
    else:    
        plt.plot(df.logmstar.loc[sel], df.EBV_excess.loc[sel],
            marker[key], ms = 10, label = labels[key], zorder = zo, alpha = 0.5, mew = 0)

xmin = 8#np.nanmin(jhu.logmstar)
xmax = 11.75#np.nanmax(jhu.logmstar)
dy = 0.25
bins = np.arange(xmin,xmax,dy)
bin_center = (bins[0:-1]+bins[1:])/2
ebv_50p = np.zeros(len(bins)-1)
ebv_16p = np.zeros(len(bins)-1)
ebv_84p = np.zeros(len(bins)-1)
ebv_med = np.zeros(len(bins)-1)
for i in range(len(bins)-1):
    subset = df.EBV_excess[(np.logical_and(jhu.logmstar>=bins[i], jhu.logmstar<=bins[i+1]))]
    #print bins[i], len(subset)
    if len(subset>0):
        ebv_med[i] = np.nanmedian(subset)
        ebv_50p[i] = np.percentile(subset,50)
        ebv_16p[i] = np.percentile(subset,16)
        ebv_84p[i] = np.percentile(subset,84)
    else:
        ebv_med[i] = np.nan
        ebv_50p[i] = np.nan
        ebv_16p[i] = np.nan
        ebv_84p[i] = np.nan

plt.plot(bin_center,ebv_50p,'o-',color = 'orange', zorder = 5, lw = 2)
plt.plot(bin_center,ebv_16p,'o-',color = 'orange',zorder = 5, lw = 2)
plt.plot(bin_center,ebv_84p,'o-',color = 'orange', zorder = 5, lw = 2)
y = np.nanmedian(df.EBV_excess)
plt.plot((xmin,xmax),(y,y), color = 'k', lw = 3)
y=np.nanmedian(df.EBV_excess[jhu.logmstar < 9.5])
plt.plot((xmin,9.5),(y,y),  color = 'k', ls = '--', lw = 3)
y=np.nanmedian(df.EBV_excess[jhu.logmstar > 9.5])
plt.plot((9.5,xmax),(y,y), color = 'k', ls = '--', lw = 3)
plt.axvline(x=9.5)
#plt.plot(bin_center,ebv_med,'o-',color = 'red', zorder = 5, lw = 2)
plt.xlabel('Stellar Mass')
plt.ylabel('Balmer decrement [ (Halpha flux/Hbeta flux)/2.86 ]')
