# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 11:17:19 2019

@author: mugdhapolimera
"""

import pandas as pd
import numpy as np

res = pd.read_csv(r'../RESOLVE_full_blend_dext_new.csv')
res.index = res.name
res = res[['econame']]
veronagn = list(pd.read_csv(r'../../xray/catalog_matching/VeronAgnMatched.csv')['eco+res_name'])
hmqagn = list(pd.read_csv(r'../../xray/catalog_matching/HMQAgnMatched.csv')['eco+res_name'])
broadagn = list(pd.read_csv(r'../../xray/catalog_matching/BroadlineAgnMatched.csv')['eco+res_name'])
xray = list(pd.read_csv(r'../../xray/ECO+RESOLVE_xray_new.csv')['name'])
lowmass = list(pd.read_csv(r'../../xray/catalog_matching/LowMassAGNMatched.csv')['eco+res_name'])

#nuggetAGN = ['rf0013', 'rf0006', 'rf0073', 'rf0503', 'rf0127', 'rs0022', 
#             'rs1036', 'rf0190', 'rf0002']

nuggetAGN = list(pd.read_csv('carrnugget.csv')['name'])
#nuggetAGN = list(pd.read_csv('p19bluenugget.csv')['name'])

heiiagn = ['rf0013','rf0372', 'rs0006','rs0463','rs1214']
print ('Veron AGN: ',[x for x in nuggetAGN if (x in veronagn) or \
                      (res.loc[x]['econame'] in veronagn)])
print ('HMQ AGN: ',[x for x in nuggetAGN if (x in hmqagn) or \
                    (res.loc[x]['econame'] in hmqagn)])
print ('Broadline AGN: ',[x for x in nuggetAGN if (x in broadagn) or \
                          (res.loc[x]['econame'] in broadagn)])
print ('X-ray AGN: ',[x for x in nuggetAGN if (x in xray) or \
                      (res.loc[x]['econame'] in xray)])
print ('Low Mass AGN: ',[x for x in nuggetAGN if (x in lowmass) or \
                      (res.loc[x]['econame'] in lowmass)])
print ('HeII AGN: ',[x for x in nuggetAGN if (x in heiiagn) or \
                      (res.loc[x]['econame'] in heiiagn)])