# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 06:23:07 2020

@author: mugdhapolimera

Create confidence levels for SFing-AGN
"""
import numpy as np
import pandas as pd

master = pd.read_csv('ECO+RESOLVE_snr5_master.csv')
master.index = master.name
flags = pd.read_csv('eco+resolve_emlineclass_snr5_master.csv')
flags.index = flags.galname
num = len(master)
conf = pd.DataFrame({'name': master.index.values, 'nsa': np.zeros(num),
                           'jhu': np.zeros(num), 'port': np.zeros(num), 
                           'confidence_level': np.zeros(num)})

nsaagn = master.name[(master.source == 'nsa') & ((flags.sftoagn) | 
        (flags.composite) | (flags.defagn))]
jhuagn = master.name[(master.source == 'jhu') & ((flags.sftoagn) | 
        (flags.composite) | (flags.defagn))]
portagn = master.name[(master.source == 'port') & ((flags.sftoagn) | 
        (flags.composite) | (flags.defagn))]
nsasf = master.name[(master.source == 'nsa') & ((flags.defstarform))]
jhusf = master.name[(master.source == 'jhu') & ((flags.defstarform))]
portsf = master.name[(master.source == 'port') & ((flags.defstarform))]

conf.index = conf.name
conf['nsa'].loc[nsaagn] = 1
conf['nsa'].loc[nsasf] = -1

conf['port'].loc[portagn] = 1
conf['port'].loc[portsf] = -1

conf['jhu'].loc[jhuagn] = 1
conf['jhu'].loc[jhusf] = -1

conf['confidence_level'] = conf['nsa'] + conf['port'] + conf['jhu']

final_sample = master.index.values

resconf = pd.read_csv('Confidence_RESOLVE.csv')
ecoconf = pd.read_csv('Confidence_ECO.csv')
fullconf = resconf.append(ecoconf)
fullconf.index = fullconf.name
finalconf = fullconf.loc[master.name]

#finalconf = conf

#print(len(np.where(finalconf.confidence_level == 1.0)[0]), 
#      len(np.where(finalconf.confidence_level == 2.0)[0]), 
#      len(np.where(finalconf.confidence_level == 3.0)[0]))

dwarf = master.logmstar < 9.5
agn = np.unique(list(jhuagn) + list(nsaagn) + 
                list(portagn))
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

#finalconf.to_csv('Confidence_ECO+RESOLVE.csv')
#finalsfagn = finalconf.loc[unique]
#sfagn_targets = list(finalsfagn.iloc[np.where(finalsfagn['confidence_level'] > 0)].name)
#targetlist = master.loc[sfagn_targets][['radeg','dedeg']]
#sami = ['rs0010', 'rs0775']
#targetlist = targetlist.drop(sami)

