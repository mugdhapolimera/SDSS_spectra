#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 10:14:22 2019

@author: mugpol

Script to search cielo for galaxies that have been reduced by SOAR

The list of folders with reduced data is in the folder 
/srv/one/resolve/analysis_code/emlineanalysis/broademline/

...and given by the file
newbroad_dirs.batch1.txt

The list of folders with reduced data that might be problematic is given at
newbroad_dirs.txt

This script is looking for the SFing-AGN galaxies that are classified by the 
new emission line classification scheme. The SFing-AGN galaxies have a 'True'
flag next to their names in the file 'resolve_emlineclass_new.csv'

"""

import numpy as np
import os
import pandas as pd
import re


def reducedsoargals(gals):

    #folder which has the list of reduced SOAR folders 
    folder = '/srv/one/resolve/analysis_code/emlineanalysis/broademline/'
    os.chdir(folder)

    #File with the list of folder names of reduced SOAR galaxies
    reduced = np.genfromtxt('newbroad_dirs.batch1.txt', dtype = [('folder','S100'),
                                                                ('type','S1')])
    #File with the list of folder names of reduced SOAR galaxies, but are problematic
    problem = np.genfromtxt('newbroad_dirs.txt', dtype = [('folder','S100'),
                                                                ('type','S1')])

    reduced_gals = []
    for folders in reduced['folder']:
        for root, dirs, files in os.walk(folders):
            for galname in gals:
                for filename in files:
                    if re.match('linslo\S+'+galname+'gspec.fits', filename):
                        reduced_gals.append((galname,
                                             os.path.join(root, filename)))
    return np.array(reduced_gals)

#List of galaxy names to search for
flags = pd.read_csv('~/github/SDSS_spectra/resolve_emlineclass_filter_new.csv')
sfagn = list(flags.galname[flags.sftoagn])

reduced_gals = reducedsoargals(sfagn)

print(reduced_gals)
#np.savetxt('reduced_sfagn_soar.txt', reduced_gals, fmt = '%s')