# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 14:54:24 2019

@author: mugdhapolimera

Mid-IR AGN Selection using WISE IR colours

as prescribed by Sartori et al. 2015
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

os.chdir('C:\Users\mugdhapolimera\github\SDSS_spectra')
data = pd.read_csv('ECO+RESOLVE_filter_new.csv')

data