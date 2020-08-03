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
eco = 0
resolve = 1
sdsscat = 'jhu'
saveoutput = 1
hasnr = 0
bpt1snr = 0
cut = 5
if eco: 
    if sdsscat =='nsa':
        inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/NSA_ECO_crossmatched.csv'
    else:
        inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_full_raw.csv'
    #NSA_ECO_full.csv'#
    #outputfile = "C:/Users/mugdhapolimera/github/SDSS_Spectra/ECO_full_snr5_nsa.csv"
    catname = 'ECO'
if resolve:
    inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_raw.csv'
    #if portsmouth: 
        #outputfile = "C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_snr5_port.csv"
    #if jhu: 
        #outputfile = "C:/Users/mugdhapolimera/github/SDSS_Spectra/RESOLVE_full_hasnr5_jhu.csv"
    if sdsscat =='nsa': 
        inputfile = 'C:/Users/mugdhapolimera/github/SDSS_Spectra/NSA_RESOLVE_crossmatched.csv'#NSA_RESOLVE_full.csv'
    catname = 'RESOLVE'
outputfile = "C:/Users/mugdhapolimera/github/SDSS_Spectra/"+catname+\
                    "_full_hasnr5_"+sdsscat+".csv"

df = catfilter(eco, resolve, sdsscat, \
              saveoutput, hasnr, bpt1snr, inputfile, outputfile, cut)

inputfile = outputfile
#flag = 'jhuname'
#foremwfile = 'RESOLVE_full_snr5_port_foreground_fiber.csv'
#extvalfile = 'RESOLVE_SDSS_JHUextvals.csv'
#foredf = mw_dext(flag, inputfile, foremwfile, extvalfile)

#inputfile = foremwfile
flag = 'smc'
smcfile = catname+'_full_hasnr5_smcdext_'+sdsscat+'.csv'
smcdf = internal_dext(flag, sdsscat, inputfile, smcfile)

flag = 'mw'
mwfile = catname+'_full_hasnr5_mwdext_'+sdsscat+'.csv'
mwdf = internal_dext(flag, sdsscat, inputfile, mwfile)

outputfile = catname+'_full_hasnr5_dext_'+sdsscat+'.csv'
ext_corr = dext(smcfile, mwfile, outputfile)