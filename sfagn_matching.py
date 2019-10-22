# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 14:24:03 2019

@author: mugdhapolimera
"""

import pandas as pd
import os
os.chdir('C:\Users\mugdhapolimera\Desktop')

jhu = pd.read_csv('JHU_sfagn.txt',delimiter = ' ').index.levels[0]
port = pd.read_csv('Port_AGN.txt',delimiter = ' ').index.levels[0]
nsa = pd.read_csv('NSA_sfagn.txt',delimiter = ' ').index.levels[0]

print "Name \t JHU \t Port \t NSA \t Total"
for x in jhu:
    print x, "\t", (x in jhu), "\t", (x in port), "\t", (x in nsa),"\t", ((x in nsa) \
                    + (x in jhu) + (x in port))
