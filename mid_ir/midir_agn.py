# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 14:54:24 2019

@author: mugdhapolimera

Mid-IR AGN Selection using WISE IR colours

as prescribed by Sartori et al. 2015 using the method of Jarrett et al 2015, 
Stern et al 2012; and the method of Satyapal et al. 2018 
"""

import numpy as np
import pandas as pd
pd.set_option('display.max_rows',500)
import matplotlib.pyplot as plt
import os
from scipy.io.idl import readsav
import sys
#os.chdir('C:\Users\mugdhapolimera\github\SDSS_spectra')
def stern(x):
    return 0.8*np.ones(len(x))

def jarretx(y):
    return [2.2*np.ones(len(y)), 4.2*np.ones(len(y))]

def jarrety(x):
    return [1.7*np.ones(len(x)), 0.1*x+0.38]

def satyapalx(x):
    return 0.52*np.ones(len(x))

def satyapaly(x):
    return 5.78*x -24.50
if sys.platform == 'linux2':
    os.chdir('/afs/cas.unc.edu/users/m/u/mugpol/github')
else:
    os.chdir('C:/Users/mugdhapolimera/github')

##############################################################################
#Reading in RESOLVE catalogs and new WISE photometry and setting up the data
##############################################################################
path = os.path.join(os.getcwd(),'SDSS_spectra\/')
resolve = pd.read_csv(path+'RESOLVE_full_blend_dext_new.csv')
resolve.index = resolve.name

path = os.path.join(os.getcwd(),'izi\/')
#Reading in the RESOLVE internal database files
catalog = readsav(path+'resolvecatalogphot.dat')
catalog2 = readsav(path+'resolvecatalog.dat')

#Create a dictionary of mid-IR magnitudes, K-band magnitudes, Surface brightness
d = {'name' : catalog2['name'],
              'mw1' : catalog['mw1'],
              'mw2' : catalog['mw2'],
              'mw3' : catalog['mw3'],
              'mw4' : catalog['mw4'],
              'emw1' : catalog['emw1'],
              'emw2' : catalog['emw2'],
              'emw3' : catalog['emw3'],
              'emw4' : catalog['emw4'],
              'kmag' : catalog['kmag'],
              'ukidsskmag' : catalog['ukidsskmag'],
              'ukidsskflag' : catalog['ukidsskflag'],
              'mur50' : catalog['mur50']}
#Changing byter order from IDL dtype format to pandas default dtype
for x in d.keys():
    d[x] = d[x].byteswap().newbyteorder() 
#Making a pandas dataframe using the dictionary
df = pd.DataFrame(data = d)
df.index = df.name

#Reading in the new WISE photometry values and converting into a pandas DF
path = os.path.join(os.getcwd(),'SDSS_spectra\/')
wise = readsav(path+'resolve_wise_102919.dat')
wisedf = pd.DataFrame.from_records(data = wise['resolve_wise'][0][0])
wisedf = wisedf.apply(lambda x: x.astype(str(x.dtype).replace('>','')))
wisedf.columns = [x.lower() for x in wisedf.columns.values]
wisedf.index = wisedf.name

#Update the original DF with the new photometry
#SKIP if you want to use the original photometry
df.update(wisedf)

##############################################################################
#Performing quality control on the data
##############################################################################

#Removing nans, 0 and applying S/N > 5 thresholding for the mid-IR mags
baderr = np.isnan(df.emw1) | np.isnan(df.emw2) | np.isnan(df.emw3) | \
        np.isnan(df.emw4)
badphot = (df.mw1 == 0.0) & (df.mw2 == 0.0) & (df.mw3 == 0.0) & \
            (df.mw4 == 0.0)
snr = (df.mw1/df.emw1 > 5) & (df.mw2/df.emw2 > 5) & \
        (df.mw3/df.emw3 > 5) & (df.mw4/df.emw4 > 5)
bad = baderr | badphot | ~snr
df = df[~bad] #Removing bad data from the DF

#Checking UKIDSS and 2MASS k-band magnitudes- 
#both k-band mags > 0 
df.ukidsskflag = [int(x) for x in df.ukidsskflag]
good_kmag = (~df['ukidsskflag'] & (df['ukidsskmag'] > 0.0)) & (df['kmag'] > 0) 

#UKIDSS and 2MASS mags should be within 10 percent of each other 
#since they have wavelength bands
kmag_agree = (df['ukidsskmag']/df['kmag'] > 0.90) & \
            (df['ukidsskmag']/df['kmag'] < 1.1)

kmagflags = good_kmag & kmag_agree
df = df[kmagflags] #Removing data with bad kmags from DF

##############################################################################
#AGN classification based on WISE colour magnitudes
##############################################################################
w12 = df['mw1'] - df['mw2']
w23 = df['mw2'] - df['mw3']
w12_err = np.sqrt(df['emw1']**2 + df['emw2']**2)
w23_err = np.sqrt(df['emw2']**2 + df['emw3']**2)
big_err = (w23_err > 0.5) | (w12_err > 0.5)
big_err_name = df.name[big_err]

df = df[~big_err]
w12 = w12
w23 = w23
w12_err = w12_err
w23_err = w23_err
#mid-IR AGN if data+/-error satisfies the AGN criteria
#Stern OR Jarrett
#midiragn = ((w12 >= 0.8) | ((w23 > 2.2) & (w23 < 4.2) & (w12 < 1.7) & 
#                            (0.1*w23 + 0.38 < w12)) | 
#           (w12 >= 0.52) & (w12 >= (5.78*w23) -24.50))

#Stern OR Jarrett OR Satyapal
midiragn = ((w12-w12_err >= 0.8) | ((w23 > 2.2) & (w23 < 4.2) & \
             (w12-w12_err < 1.7) & (0.1*w23 + 0.38 < w12-w12_err)) | \
           (w12-w12_err >= 0.52) & (w12-w12_err >= (5.78*w23) -24.50))
#Stern AND Jarret AND Satyapal
#midiragn = (((w12-w12_err) >= 0.8) & (((w23-w23_err) > 2.2) & ((w23+w23_err) < 4.2) & \
#             (w12-w12_err < 1.7) & (0.1*(w23-w23_err) + 0.38 < w12-w12_err)) & \
#           ((w12-w12_err) >= 0.52) & ((w12-w12_err) >= (5.78*w23) -24.50))

#plt.plot(w23[midiragn],w12[midiragn],'rs')
#plt.ylim(-6.0,10)

plt.figure()
plt.plot(w23,w12,'o', label = 'Galaxies with reliable WISE mags')
xaxis = np.linspace(min(w23),max(w23))
#yaxis = np.linspace(min(w12), max(w12))
yaxis = np.linspace(jarrety(np.array([2.2]))[1],1.7)
plt.plot(xaxis, stern(xaxis), 'k-.', label = 'Stern12')
plt.plot(xaxis, satyapalx(xaxis), 'k', label = 'Satyapal')
xaxis = np.linspace(4.3,max(w23))
plt.plot(xaxis, satyapaly(xaxis), 'k')

xaxis = np.linspace(2.2,4.2)
plt.plot(jarretx(yaxis)[0], yaxis, 'k--', jarretx(yaxis)[1], yaxis, 'k--')
plt.plot(xaxis, jarrety(xaxis)[0], 'k--')
plt.plot(xaxis, jarrety(xaxis)[1],'k--', label = 'Jarrett15')
plt.xlabel('W2 - W3')
plt.ylabel('W1 - W2')
plt.ylim(min(w12)-0.1, max(w12)+0.1)

plt.errorbar(w23[midiragn],w12[midiragn],fmt = 'rs', xerr = w23_err[midiragn],
             yerr = w12_err[midiragn], label = 'Mid-IR AGN')
plt.legend()

midiragnname = df.name[midiragn]
print(resolve.loc[midiragnname][['radeg','dedeg']])    
print('{} mid-IR AGN out of {} galaxies having reliable WISE mags : {}%'\
      .format(len(midiragnname), len(df), \
              round(len(midiragnname)*100.0/len(df),2)))
#Flags to check SFing-AGN    
#flags = pd.read_csv('SDSS_spectra/resolve_emlineclass_full_snr5.csv')
#flags.index = flags.galname
#sfagn = list(flags.galname.iloc[np.where(flags.sftoagn)])
#sfagnndx = [x for x in range(len(resname)) if resname[x] in sfagn]
#plt.plot(w23[sfagnndx],w12[sfagnndx],'gs', label = 'SFing-AGN')

#sfagnmidir = [resname[midiragn][x] for x in range(len(resname[midiragn])) if resname[midiragn][x] in sfagn]
#print sfagnmidir
