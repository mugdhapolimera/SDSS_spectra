# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 12:46:19 2020

@author: mugdhapolimera

SDSS Catalog filtering and De-extinction Wrapper
"""

import numpy as np
import pandas as pd

#dext = __import__("SMC+MW_dext")
from catalog_filter import catfilter
from SMCandMW_dext import dext
from extinction_correction import internal_dext
from extinction_correction import mwforeground_dext as mw_dext

eco = 1
resolve = 0
full = 0
sdsscat = 'port'
saveoutput = 1
hasnr = 0
bpt1snr = 0
he2 = 0
cut = 5
he2cut = 0

if eco: 
    if sdsscat =='nsa':
        inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/NSA_ECO_crossmatched.csv'
    else:
        inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_full_raw.csv'
    survey = 'ECO'

if resolve:
    inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_raw.csv'

    if sdsscat =='nsa': 
        inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/NSA_RESOLVE_crossmatched.csv'
    survey = 'RESOLVE'

if full: #DO NOT USE FULL OPTION FOR sdsscat='port'. ONE FILE DOES NOT EXIST FOR ECO+RESOLVE
        inputfile = 'ECO+RESOLVE_full_raw.csv'
        survey = 'ECO+RESOLVE'

#Modify the outputfile name according to how you filter the sample
outputfile = "C:/Users/mugdhapolimera/github/SDSS_Spectra/"+survey+\
                    "_full_snr5"+sdsscat+".csv"
df = catfilter(eco, resolve, sdsscat, \
              saveoutput, hasnr, bpt1snr, he2, inputfile, outputfile, cut, he2cut)

inputfile = outputfile

if sdsscat == 'port':
    flag = 'portname'
    foremwfile = survey+'_full_snr5_'+sdsscat+'_foreground_fiber.csv'
    extvalfile = survey+'_SDSS_JHUextvals.csv'
    foredf = mw_dext(flag, inputfile, foremwfile, extvalfile)

    inputfile = foremwfile


flag = 'smc'
smcfile = survey+'_full_snr5_SMCdext_'+sdsscat+'.csv'
smcdf = internal_dext(flag, sdsscat, inputfile, smcfile)

flag = 'mw'
mwfile = survey+'_full_snr5_MWdext_'+sdsscat+'.csv'
mwdf = internal_dext(flag, sdsscat, inputfile, mwfile)


#Final output file with 
outputfile = survey+'_full_snr5_dext_'+sdsscat+'.csv'
ext_corr = dext(smcfile, mwfile, outputfile)
